function figure_qbold_nonoise(simdir)

	rng('default'); % reset random number generator

	TEGRE=[34e-3 36e-3]';
	TEGREi=35e-3;
	TE=80e-3;
	tauASE=[0 (12:4:68)]'./1000;
	tauASEi=[0 (15:3:66)]'./1000;
	tau_cutoff=15e-3;

	Ds=[120 60 30 20 10 5.6 15 30 45 90 180];
	Rs=Ds./2;
	
	aVessels=pi.*Rs.^2;
	lVessels=[5390 2690 1350 900 450 600 450 900 1350 2690 5390];
	nVessels=[1880 1.5e4 1.15e5 3.92e5 3.01e6 5.92e7 3.01e6 3.92e5 1.15e5 1.5e4 1880];

	volVessels=nVessels.*lVessels.*aVessels;
	relVf=volVessels./sum(volVessels);
	
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf3pc_sharan/simvessim_res' num2str(Rs(k)) '.mat']);
		sppR(:,:,k)=spp;
		pR(k)=p;
	end	
	
	N=1000;
	Vtot=rand(N,1).*0.1; %Total CBV between 1 and 10%
	
	Vf=repmat(relVf,N,1).*repmat(Vtot,1,length(Rs));
	
	PaO2=rand(N,1).*(130-100)+100; %Baseline PaO2 between 100 and 130mmHg
	dPaO2=rand(N,1).*(320-270)+270; %Delta PaO2 between 270 and 320mmHg
	%PaO2=ones(N,1).*115;
	%dPaO2=ones(N,1).*295; %Delta PaO2 between 270 and 320mmHg
	%PaO2=[ones(N/2,1).*130; ones(N/2,1).*100]; %Baseline PaO2 between 100 and 130mmHg
	%dPaO2=[ones(N/2,1).*320; ones(N/2,1).*270]; %Delta PaO2 between 270 and 320mmHg
	%dPaO2=[ones(N/4,1).*370; ones(N/4,1).*320; ones(N/4,1).*270; ones(N/4,1).*220];
	PaO2h=PaO2+dPaO2;
	
	Ya=1./(23400./(PaO2.^3+150.*PaO2)+1);
	Yah=1./(23400./(PaO2h.^3+150.*PaO2h)+1);
	E0=rand(N,1); %Baseline OEF between 0 and 100%
	
	Hct=(rand(N,1).*(0.5-0.37)+0.37).*0.85;
	Hct_assumed=0.34;
	
	phi=1.34;
	epsilon=0.0031;
	Hb=Hct./0.03;

	CaO2=phi.*Hb.*Ya+PaO2.*epsilon;
	CaO2h=phi.*Hb.*Yah+PaO2h.*epsilon; 

	CmetO2=CaO2.*E0; 

	PvO2=calcPvO2(CaO2-CmetO2,Hb); 
	PvO2h=calcPvO2(CaO2h-CmetO2,Hb);

	Yv=1./(23400./(PvO2.^3+150*PvO2)+1);
	Yvh=1./(23400./(PvO2h.^3+150*PvO2h)+1);

	k=0.4;
	Yc=Ya.*k+Yv.*(1-k);
	Ych=Yah.*k+Yvh.*(1-k);

	Y=[Ya Ya Ya Ya Ya Yc Yv Yv Yv Yv Yv];
	Yh=[Yah Yah Yah Yah Yah Ych Yvh Yvh Yvh Yvh Yvh];
	
	for j=1:N
		for k=1:length(Rs)
			[sigASE(:,k) tauASE sigASEev(:,k) sigASEiv(:,k)]=generate_signal(pR(k),sppR(:,:,k),'display',false,'Vf',Vf(j,k),'Y',Y(j,k),'Hct',Hct(k),'seq','ASE','includeIV',true,'T2EV',80e-3,'T2b0',189e-3,'TE',TE,'tau',tauASE);
			sigASEev(:,k)=sigASEev(:,k)./exp(-TE./80e-3);
			
			[sigGREBASE(:,k) teGRE sigGREevBASE(:,k) sigGREivBASE(:,k)]=generate_signal(pR(k),sppR(:,:,k),'display',false,'Vf',Vf(j,k),'Y',Y(j,k),'Hct',Hct(k),'seq','GRE','includeIV',true,'T2EV',80e-3,'T2b0',189e-3,'TE',TEGRE);
			sigGREevBASE(:,k)=sigGREevBASE(:,k)./exp(-TEGRE./80e-3);
			[sigGREHYPE(:,k) teGRE sigGREevHYPE(:,k) sigGREivHYPE(:,k)]=generate_signal(pR(k),sppR(:,:,k),'display',false,'Vf',Vf(j,k),'Y',Yh(j,k),'Hct',Hct(k),'seq','GRE','includeIV',true,'T2EV',80e-3,'T2b0',189e-3,'TE',TEGRE);
			sigGREevHYPE(:,k)=sigGREevHYPE(:,k)./exp(-TEGRE./80e-3);
		end
		
		sigASEtot(:,j)=(1-sum(Vf(j,:),2)).*prod(sigASEev,2).*exp(-TE./80e-3)+sum(bsxfun(@times,Vf(j,:),sigASEiv),2);
		sigASEtoti(:,j)=interp1(tauASE,sigASEtot(:,j),tauASEi);

		sigGREtotBASE(:,j)=(1-sum(Vf(j,:),2)).*prod(sigGREevBASE,2).*exp(-TEGRE./80e-3)+sum(bsxfun(@times,Vf(j,:),sigGREivBASE),2);
		sigGREtotBASEi(:,j)=interp1(TEGRE,sigGREtotBASE(:,j),TEGREi);
		sigGREtotHYPE(:,j)=(1-sum(Vf(j,:),2)).*prod(sigGREevHYPE,2).*exp(-TEGRE./80e-3)+sum(bsxfun(@times,Vf(j,:),sigGREivHYPE),2);
		sigGREtotHYPEi(:,j)=interp1(TEGRE,sigGREtotHYPE(:,j),TEGREi);
		
		[paramsASEtot(:,j) paramsASEtotsd(:,j)]=calc_qbold_params(p,sigASEtoti(:,j),tauASEi,tau_cutoff,Hct_assumed);
		%fprintf([num2str(j) '.']);
		disp(num2str(j));
	end
		
	TEGRE=TEGREi;
	const=4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*Hct; %qBOLD model scaling constants
	scale=(27e-3./TEGRE+0.2).*(245.1./dPaO2+0.1); %Scale factor from Blockley et al., 2013
	hqDBV=(sigGREtotHYPEi'./sigGREtotBASEi'-1).*scale;
	hqOEF=paramsASEtot(1,:)'./hqDBV./const;
	sqDBV=paramsASEtot(2,:)';
	sqOEF=paramsASEtot(3,:)';
	R2p=paramsASEtot(1,:)';
	rho=104; %density of tissue in g/dl^tissue
	dhb_content=(sum(Vf(:,6:11),2)./rho).*Hb.*E0.*100;
	
	trueDBV1=paramsASEtot(1,:)'./(const.*E0); %Based on DBV required to estimate true OEF
	trueDBV2=sum(Vf(:,7:11),2); %Based on sum of vessels with dHb (no capillary, only venules -> veins)
	trueR2p=const.*E0.*trueDBV2;
	
	%FIGURE 2A	
	figure;
	scatter(trueR2p,R2p,[],dhb_content,'filled');
	colormap parula;
	hold on;
	plot([0 14],[0 14],'k:');
	axis square;
	box on;
	grid on;
	title('Fig. 2a. Apparent R2p vs SDR R2p')
	ylim([0 14]);
	ylabel('SDR R2p (Hz)')
	xlim([0 14]);
	xlabel('Apparent R2p (Hz)')
	colorbar

    XR2p=[ones(size(trueR2p)) trueR2p];
    aR2p=lscov(XR2p,R2p);
    disp(['R2p gradient = ' num2str(round(aR2p(2),2)) ', intercept = ' num2str(round(aR2p(1),2))]);

	%FIGURE 2B	
	figure;
	scatter(trueDBV2.*100,sqDBV.*100,[],E0.*100,'filled');
	colormap parula;
	hold on;
	plot([-1 5],[-1 5],'k-');
	axis square;
	box on;
	grid on;
    set(gca,'XTick',[0:1:5])
    set(gca,'YTick',[0:1:5])
	title('Fig. 2b. sqBOLD DBV')
	ylim([-0.5 5]);
	ylabel('True DBV (%)')
	xlim([-0.5 5]);
	xlabel('Apparent DBV (%)')	
	colorbar
	
	%FIGURE 2C	
	figure;
	scatter(trueDBV2.*100,hqDBV.*100,[],E0.*100,'filled');
	colormap parula;
	hold on;
	plot([-1 5],[-1 5],'k-');
	axis square;
	box on;
	grid on;
    set(gca,'XTick',[0:1:5])
    set(gca,'YTick',[0:1:5])
	title('Fig. 2c. hqBOLD DBV')
	ylim([-0.5 5]);
	ylabel('True DBV (%)')
	xlim([-0.5 5]);
	xlabel('Apparent DBV (%)')	
	colorbar

    XhqDBV=[ones(size(trueDBV2)) trueDBV2.*100];
    ahqDBV=lscov(XhqDBV,hqDBV.*100);
    %restricted range (rr) where OEF>30%
    XhqDBVrr=[ones(size(trueDBV2(E0>0.3))) trueDBV2(E0>0.3).*100];
    ahqDBVrr=lscov(XhqDBVrr,hqDBV(E0>0.3).*100);
    disp(['hqDBV gradient = ' num2str(round(ahqDBV(2),2)) ', intercept = ' num2str(round(ahqDBV(1),2))]);
    disp(['hqDBV (OEF>30%) gradient = ' num2str(round(ahqDBVrr(2),2)) ', intercept = ' num2str(round(ahqDBVrr(1),2))]);

	%FIGURE 2D
	figure;
	cmap=colormap('lines');
	scatter(E0.*100,(sqDBV-trueDBV2)./trueDBV2.*100,[],cmap(1,:),'filled');
	hold on;
	scatter(E0.*100,(hqDBV-trueDBV2)./trueDBV2.*100,[],cmap(2,:),'filled');
	legend('sqBOLD DBV','hqBOLD DBV','location','northwest')
	axis square;
	box on;
	grid on;
    set(gca,'XTick',[0:20:100])
    set(gca,'YTick',[-100:100:300])
	title('Fig. 2d. Error in apparent DBV')
    ylim([-150 350])
	ylabel('Error in apparent DBV (%)')
    xlim([0 100])
	xlabel('True OEF (%)')	

	%FIGURE 2E
	figure;
	scatter(E0.*100,sqOEF.*100,[],trueDBV2.*100,'filled');
	colormap parula;
	hold on;
	plot([0 100],[0 100],'k-');
	axis square;
	box on;
	grid on;
    set(gca,'XTick',[0:20:100])
    set(gca,'YTick',[0:20:100])
	title('Fig. 2e. sqBOLD OEF')
	ylim([0 100]);
	ylabel('True OEF (%)')
	xlim([0 100]);
	xlabel('Apparent OEF (%)')
	colorbar
	
	%FIGURE 2F
	figure;
	scatter(E0.*100,hqOEF.*100,[],trueDBV2.*100,'filled');
	colormap parula;
	hold on;
	plot([0 100],[0 100],'k-');
	axis square;
	box on;
	grid on;
    set(gca,'XTick',[0:20:100])
    set(gca,'YTick',[0:20:100])
	title('Fig. 2f. hqBOLD OEF')
	ylim([0 100]);
	ylabel('True OEF (%)')
	xlim([0 100]);
	xlabel('Apparent OEF (%)')
	colorbar
	
	function PvO2=calcPvO2(CvO2,Hb)
	%return venous partial pressure of oxygen (PvO2) based on input of venous oxygen content
	%and haemoglobin concentration
	
	CvO2=CvO2(:);
	N=length(CvO2);

	phi=1.34;
	epsilon=0.0031;
	%a=phi.*Hb;
	b=epsilon;
	
	for k=1:N
	
		if CvO2(k)==0
			PvO2(k,:)=0;
			return;
		end

		a=phi.*Hb(k);

		A=b;
		CvO2temp=CvO2(k);
		B=a-CvO2temp;
		C=150*b;
		D=150*a+23400*b-150*CvO2temp;
		E=-23400*CvO2temp;

		rts=roots([A B C D E]);
		PvO2temp=rts((rts>0) & (imag(rts)==0));

		if isempty(PvO2temp)
			PvO2(k,:)=NaN;
		else
			PvO2(k,:)=PvO2temp;
		end
		
	end

	return