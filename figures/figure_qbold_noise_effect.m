function figure_qbold_noise_effect(simdir)

	rng('default'); % reset random number generator

	TEGRE=36e-3;
	TE=80e-3;
	tauASE=[0 (12:4:68)]./1000;
	%tauASEi=[0 (16:4:64)]./1000;
	tauASEi=[0 (15:3:66)]./1000;
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
	
	N=20000;
	Vtot=ones(N,1).*0.05; %Fixed value of total CBV at 5%
	
	Vf=repmat(relVf,N,1).*repmat(Vtot,1,length(Rs));
	
	PaO2=ones(N,1).*110; %Fixed value of baseline PaO2 at 110mmHg
	dPaO2=ones(N,1).*290; %Fixed value of change in PaO2 at 290mmHg
	PaO2h=PaO2+dPaO2;
	
	Ya=1./(23400./(PaO2.^3+150.*PaO2)+1);
	Yah=1./(23400./(PaO2h.^3+150.*PaO2h)+1);
	E0=ones(N,1).*0.4; %Fixed value of OEF at 40%
		
	Hct=0.34;
	
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
	
	ASESNR=52; %Image SNR
	BOLDSNR=87;
	nBOLD=100; %Number of baseline BOLD measurements
	nBOLDh=100; %Number of hyperoxia BOLD measurements

	for k=1:length(Rs)
		[sigASE(:,k) tauASE sigASEev(:,k) sigASEiv(:,k)]=generate_signal(pR(k),sppR(:,:,k),'display',false,'Vf',Vf(1,k),'Y',Y(1,k),'Hct',Hct,'seq','ASE','includeIV',true,'T2EV',80e-3,'T2b0',189e-3,'TE',TE,'tau',tauASE);
		sigASEev(:,k)=sigASEev(:,k)./exp(-TE./80e-3);
		
		[sigGREBASE(:,k) teGRE sigGREevBASE(:,k) sigGREivBASE(:,k)]=generate_signal(pR(k),sppR(:,:,k),'display',false,'Vf',Vf(1,k),'Y',Y(1,k),'Hct',Hct,'seq','GRE','includeIV',true,'T2EV',80e-3,'T2b0',189e-3,'TE',TEGRE);
		sigGREevBASE(:,k)=sigGREevBASE(:,k)./exp(-TEGRE./80e-3);
		[sigGREHYPE(:,k) teGRE sigGREevHYPE(:,k) sigGREivHYPE(:,k)]=generate_signal(pR(k),sppR(:,:,k),'display',false,'Vf',Vf(1,k),'Y',Yh(1,k),'Hct',Hct,'seq','GRE','includeIV',true,'T2EV',80e-3,'T2b0',189e-3,'TE',TEGRE);
		sigGREevHYPE(:,k)=sigGREevHYPE(:,k)./exp(-TEGRE./80e-3);
	end
			
	for j=1:N

		sigASEtot(:,j)=(1-sum(Vf(j,:),2)).*prod(sigASEev,2).*exp(-TE./80e-3)+sum(bsxfun(@times,Vf(j,:),sigASEiv),2);
		sigASEtoti(:,j)=interp1(tauASE,sigASEtot(:,j),tauASEi);
		sigASEtoti(:,j)=sigASEtoti(:,j)+sigASEtoti(1,j).*randn(size(sigASEtoti(:,j)))./ASESNR;
		
		sigGREtotBASE(:,j)=(1-sum(Vf(j,:),2)).*prod(sigGREevBASE,2).*exp(-TEGRE./80e-3)+sum(bsxfun(@times,Vf(j,:),sigGREivBASE),2);
		sigGREtotBASE(:,j)=mean(sigGREtotBASE(:,j)+sigGREtotBASE(1,j).*randn(nBOLD,1)./BOLDSNR);
		sigGREtotHYPE(:,j)=(1-sum(Vf(j,:),2)).*prod(sigGREevHYPE,2).*exp(-TEGRE./80e-3)+sum(bsxfun(@times,Vf(j,:),sigGREivHYPE),2);
		sigGREtotHYPE(:,j)=mean(sigGREtotHYPE(:,j)+sigGREtotBASE(1,j).*randn(nBOLDh,1)./BOLDSNR);
		
		[paramsASEtot(:,j) paramsASEtotsd(:,j)]=calc_qbold_params(p,sigASEtoti(:,j),tauASEi,tau_cutoff);
		%fprintf([num2str(j) '.']);
		disp(num2str(j));
	end
	
	const=4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*Hct; %qBOLD model scaling constants
	%dPaO2=PaO2h-PaO2
	scale=(27e-3./TEGRE+0.2)*(dPaO2./245.1+0.1); %Scale factor from Blockley et al., 2013
	hqDBV=(sigGREtotHYPE'./sigGREtotBASE'-1).*scale;
	hqOEF=paramsASEtot(1,:)'./hqDBV./const;
	sqDBV=paramsASEtot(2,:)';
	sqOEF=paramsASEtot(3,:)';
	
	trueDBV1=paramsASEtot(1,:)'./(const.*E0); %Based on DBV required to estimate true OEF
	trueDBV2=sum(Vf(:,7:11),2); %Based on sum of vessels with dHb (capillary -> veins)
	trueR2p=const.*E0.*trueDBV2;
	
	r2p=paramsASEtot(1,:)';

	figure;
	
	%FIGURE 3a	
	subplot(231)
	histogram(r2p,linspace(-2,6,64),'normalization','probability')
	title('Fig. 3a. Apparent R_2^\prime')
	ylabel('Probability')
	xlim([-2 6])
	ylim([0 0.2])
	xlabel('Apparent R_2^\prime (s^{-1})')
	set(gca,'xtick',[-2:2:6])
	axis square;
	box on;
	grid on;

	%FIGURE 3B
	subplot(232)
	histogram(sqDBV.*100,linspace(-4,12,64),'normalization','probability')
	title('Fig. 3b. sqBOLD DBV')
	ylabel('Probability')
	xlim([-4 12])
	ylim([0 0.05])
	xlabel('Apparent DBV (%)')
	set(gca,'xtick',[-4:4:12])
	set(gca,'ytick',[0:0.01:0.05])
	axis square;
	box on;
	grid on;	
	
	%FIGURE 3C
	subplot(233)
	histogram(hqDBV.*100,linspace(-4,12,64),'normalization','probability')
	title('Fig. 3c. hqBOLD DBV')
	ylabel('Probability')
	xlim([-4 12])
	ylim([0 0.4])
	xlabel('Apparent DBV (%)')
	set(gca,'xtick',[-4:4:12])
	axis square;
	box on;
	grid on;	
	
	%FIGURE 3D
	subplot(235)
	histogram(sqOEF.*100,linspace(-40,120,64),'normalization','probability')
	title('Fig. 3d. sqBOLD OEF')
	ylabel('Probability')
	xlim([-40 120])
	ylim([0 0.15])
	xlabel('Apparent OEF (%)')
	set(gca,'xtick',[-40:40:120])
	axis square;
	box on;
	grid on;
	
	%FIGURE 3E
	subplot(236)
	histogram(hqOEF.*100,linspace(-40,120,64),'normalization','probability')
	title('Fig. 3e. hqBOLD OEF')
	ylabel('Probability')
	xlim([-40 120])
	ylim([0 0.3])
	xlabel('Apparent OEF (%)')
	set(gca,'xtick',[-40:40:120])
	axis square;
	box on;
	grid on;		

	figure;

	%FIGURE S1A
	subplot(121)
	bcx=(0:0.001:0.1001)';
	bex=bcx-0.001/2;
	bcy=(0:0.01:1.01)';
	bey=bcy-0.01/2;
	H=histogram2(sqDBV,sqOEF,bex,bey,'facecolor','flat');
	sqBOLDhist=H.Values;
	imagesc(flipud(rot90(sqBOLDhist)))
	cmap=colormap('hot');
	cmap2=flipud(cmap);
	colormap(cmap2);
	set(gca,'ydir','normal')
	set(gca,'xtick',[1:10:101])
	set(gca,'xticklabel',bcx(1:10:101).*100)
	set(gca,'ytick',[1:10:101])
	set(gca,'yticklabel',bcy(1:10:101).*100)
	title('Fig. S1a. sqBOLD OEF vs sqBOLD DBV')
	xlabel('Apparent DBV (%)');
	ylabel('Apparent OEF (%)');
	xlim([1 101])
	ylim([1 101])
	axis square;
	box on;
	grid on;
	
	%FIGURE S1B
	subplot(122)
	bcx=(0:0.001:0.1001)';
	bex=bcx-0.001/2;
	bcy=(0:0.01:1.01)';
	bey=bcy-0.01/2;
	H=histogram2(hqDBV,hqOEF,bex,bey,'facecolor','flat');
	sqBOLDhist=H.Values;
	imagesc(flipud(rot90(sqBOLDhist)))
	cmap=colormap('hot');
	cmap2=flipud(cmap);
	colormap(cmap2);
	set(gca,'ydir','normal')
	set(gca,'xtick',[1:10:101])
	set(gca,'xticklabel',bcx(1:10:101).*100)
	set(gca,'ytick',[1:10:101])
	set(gca,'yticklabel',bcy(1:10:101).*100)
	title('Fig. S1b. hqBOLD OEF vs hqBOLD DBV')
	xlabel('Apparent DBV (%)');
	ylabel('Apparent OEF (%)');
	xlim([1 101])
	ylim([1 101])
	axis square;
	box on;
	grid on;
	
	%keyboard;
	
	function PvO2=calcPvO2(CvO2,Hb)
	%return venous partial pressure of oxygen (PvO2) based on input of venous oxygen content
	%and haemoglobin concentration
	
	CvO2=CvO2(:);
	N=length(CvO2);

	phi=1.34;
	epsilon=0.0031;
	a=phi*Hb;
	b=epsilon;
	
	for k=1:N
	
		if CvO2(k)==0
			PvO2(k,:)=0;
			return;
		end

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