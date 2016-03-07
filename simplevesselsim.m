function [storedProtonPhase p r]=simplevesselsim(p,seed)
    
    %calculate standard deviation of proton diffusion step
    p.stdDev = sqrt(2*p.D*p.dt);
    
    %set up random number generator
    if nargin>1
    	p.seed=seed;
    else
    	fid=fopen('/dev/urandom');
    	p.seed=fread(fid, 1, 'uint32');
    end
    rng(p.seed); %use a random seed to avoid problems when running on a cluster
    
for k=1:p.N
    
    %set up universe
    [p.universeSize, vesselOrigins, vesselNormals, protonPosit, r(k).numVessels, r(k).vesselVolFrac] = setupUniverse(p);

	%generate random walk path
	[protonPosits p.numSteps] = randomWalk(p,protonPosit);

	%calculate field at each point
	fieldAtProtonPosit=calculateField(p, protonPosits, vesselOrigins, vesselNormals, r(k).numVessels);
	
	%calculate phase at each point
	pts=round(p.deltaTE./p.dt); %pts per deltaTE
	storedProtonPhase(:,k)=sum(reshape(fieldAtProtonPosit,pts,p.numSteps/pts).*p.gamma.*p.dt,1)';

end

	%keyboard;

return;

% Set up the universe of cylindrical vessels
function [universeSize, vesselOrigins, vesselNormals, protonPosit, numVessels, vesselVolFrac] = setupUniverse(p)

    universeSize = 2*(sqrt(2*1.3e-9*p.TE)+40*p.R);%as in Boxerman
    volUniverse = (4/3)*pi*universeSize^3;
    M=1000; %max number of vessels
    
    %uniform random distribution of vessel seed points in sphere
    %vesselOrigins=(rand(M,3)-0.5).*2.*universeSize;
    %withinSphere=find(sqrt(sum(vesselOrigins.^2,2))<=universeSize);
    %vesselOrigins=vesselOrigins(withinSphere,:);
    
    %distribute some vessel origins within sphere and some on surface (50-50)
    randomNormals=randn(M,3);
    randomNormals=randomNormals./repmat(sqrt(sum(randomNormals.^2,2)),1,3);
    r=repmat(universeSize.*rand(M,1).^(1/3),1,3);
    %r=repmat(universeSize,M,3);
    r(2:2:end,:)=repmat(universeSize,M/2,3); %half of vessel origins on the surface
    vesselOrigins=r.*randomNormals;
    
    %uniform random distribution of random normals (orientations)
    %vesselNormals=randn(length(withinSphere),3);
    %vesselNormals=vesselNormals./repmat(sqrt(sum(vesselNormals.^2,2)),1,3);
    %vesselNormals=repmat([0 0 1],length(withinSphere),1);

    vesselNormals=randn(M,3);
    vesselNormals=vesselNormals./repmat(sqrt(sum(vesselNormals.^2,2)),1,3);
    
    %calculate lengths of vessels in sphere
    a=sum(vesselNormals.^2,2);
    b=2*sum(vesselOrigins.*vesselNormals,2);
    c=sum(vesselOrigins.^2,2)-universeSize.^2;
    delta=b.*b-4*a.*c;
    u1=(-b-sqrt(delta))./2./a;
    u2=(-b+sqrt(delta))./2./a;
    p1=vesselOrigins+repmat(u1,1,3).*vesselNormals;
    p2=vesselOrigins+repmat(u2,1,3).*vesselNormals;
    l=sqrt(sum((p2-p1).^2,2));
    
    %find vessel number cutoff for desired volume fraction
    volSum=(cumsum(l).*pi.*p.R.^2);
	cutOff=find(volSum<(volUniverse.*p.sphereDensity),1,'last');
    
    if cutOff==M
    	disp('Error: Increase max vessels!');
    end
    
    vesselOrigins=vesselOrigins(1:cutOff,:);
    vesselNormals=vesselNormals(1:cutOff,:);
	vesselVolFrac = volSum(cutOff)/volUniverse;
    numVessels = cutOff;
    
    protonPosit=[0 0 0];
    
    %figure;
    %plot3([p1(1:cutOff,1); p2(1:cutOff,1)],[p1(1:cutOff,2); p2(1:cutOff,2)],[p1(1:cutOff,3); p2(1:cutOff,3)],'-')
    
    %keyboard;
return;

%random walk
function [protonPosits numSteps] = randomWalk(p,protonPosit);
	numSteps=round((p.TE*2)/p.dt);
	protonPosits=p.stdDev.*randn(numSteps,3);
	protonPosits(1,:)=protonPosit;
	protonPosits=cumsum(protonPosits);
return;

%calculate magnetic field at proton location
function [totalField] = calculateField(p, protonPosits, vesselOrigins, vesselNormals, numVessels)
	
	numSteps=size(protonPosits,1);
	
	protonPosits=repmat(permute(protonPosits,[3 2 1]),numVessels,1,1);
	vesselOrigins=repmat(vesselOrigins,1,1,numSteps);
	vesselNormals=repmat(vesselNormals,1,1,numSteps);
	
	relPosits=protonPosits-vesselOrigins;
	
	%perpendicular distance from proton to vessel
	r=sqrt(sum((relPosits-repmat(dot(relPosits,vesselNormals,2),1,3,1).*vesselNormals).^2,2));

	%elevation angle between vessel and the z-axis (just do it for one time step and repmat later)
	theta=acos(dot(vesselNormals(:,:,1),repmat([0 0 1],numVessels,1,1),2));
	
	np=zeros(numVessels,3,numSteps);
	np=relPosits-repmat(dot(relPosits,vesselNormals,2),1,3,1).*vesselNormals;
	np=np./repmat(sqrt(sum(np.^2,2)),1,3,1);
	nb=cross(repmat([0 0 1],numVessels,1,numSteps),vesselNormals);
	nb=nb./repmat(sqrt(sum(nb.^2,2)),1,3,1);
	nc=cross(vesselNormals,nb);
	nc=nc./repmat(sqrt(sum(nc.^2,2)),1,3,1);

	%azimuthal angle in plane perpendicular to vessel
	%phi=zeros(numVessels,3,numSteps);
	phi=acos(dot(np,nc,2));
	
	fields=p.B0.*2.*pi.*p.deltaChi.*(p.R./r).^2.*cos(2.*phi).*sin(repmat(theta,1,1,numSteps)).^2;
	%fields(r<p.R)=p.B0.*2.*pi./3.*p.deltaChi.*(3.*cos(repmat(theta,1,1,numSteps)).^2-1);
	%keyboard;
	fields(r<p.R)=0; %really should be calculating field inside vessel

    totalField= p.B0 + sum(fields,1);
    totalField=squeeze(totalField);  
    
    %keyboard;
return;


% Calculate phase at proton position
function [phase] = calculatePhase(p, totalField, phase)
    phase = phase + p.gamma*totalField*p.dt;
return;


