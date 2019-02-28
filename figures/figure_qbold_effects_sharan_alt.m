function figure_qbold_effects_sharan_alt(simdir)

	%acquisition parameter as in An & Lin, MRM 50:708-716 (2003)
	TE=64e-3;
	tauASE=[0 (12:4:20)]./1000;
	tau_cutoff=10e-3;

	Ds=[120 60 30 20 10 5.6 15 30 45 90 180];
	Rs=Ds./2;
	
	aVessels=pi.*Rs.^2;
	lVessels=[5390 2690 1350 900 450 600 450 900 1350 2690 5390];
	nVessels=[1880 1.5e4 1.15e5 3.92e5 3.01e6 5.92e7 3.01e6 3.92e5 1.15e5 1.5e4 1880];
	volVessels=nVessels.*lVessels.*aVessels;
	relVf=volVessels./sum(volVessels);
	Vtot=0.05;
	Vf=relVf.*Vtot;
	
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf3pc_sharan/simvessim_res' num2str(Rs(k)) '.mat']);
		sppR(:,:,k)=spp;
	end	
	
	N=1000;
	Vtotrange=[0 0.1];
	Vtot=rand(N,1).*(max(Vtotrange)-min(Vtotrange))+min(Vtotrange);
	Vf=repmat(relVf,N,1).*repmat(Vtot,1,length(Rs));
	
	Ya=0.98;
	E0range=[0 1];
	E0=rand(N,1).*(max(E0range)-min(E0range))+min(E0range);
	Yv=Ya.*(1-E0);

	k=0.4;
	Yc=Ya.*k+Yv.*(1-k);
	
	Ya=repmat(Ya,N,1);
	Y=[Ya Ya Ya Ya Ya Yc Yv Yv Yv Yv Yv];
	
	for j=1:N
		for k=1:length(Rs)
			[sigASE(:,k) tauASE sigASEev(:,k) sigASEiv(:,k)]=generate_signal(p,sppR(:,:,k),'display',false,'Vf',Vf(j,k),'Y',Y(j,k),'seq','ASE','includeIV',true,'T2EV',80e-3,'T2b0',189e-3,'TE',TE,'tau',tauASE);
			sigASEev(:,k)=sigASEev(:,k)./exp(-TE./80e-3);
		end
		sigASEtot(:,j)=(1-sum(Vf(j,:),2)).*prod(sigASEev,2).*exp(-TE./80e-3)+sum(bsxfun(@times,Vf(j,:),sigASEiv),2);
		sigASEtotCV(:,j)=(1-sum(Vf(j,7:11),2)).*prod(sigASEev(:,7:11),2).*exp(-TE./80e-3)+sum(bsxfun(@times,Vf(j,7:11),sigASEiv(:,7:11)),2);
		
		[paramsASEtot(:,j) paramsASEtotsd(:,j)]=calc_qbold_params(p,sigASEtot(:,j),tauASE,tau_cutoff);
		[paramsASEtotCV(:,j) paramsASEtotCVsd(:,j)]=calc_qbold_params(p,sigASEtotCV(:,j),tauASE,tau_cutoff);
		fprintf([num2str(j) '.']);
	end
	
	%FIGURE S6A	
	figure;
	const=4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct;
	rho=104; %density of tissue in g/dl^tissue
	hb=p.Hct./0.03; %haemoglobin concentration in g^Hb/dl^blood
	cbv=sum(Vf(:,6:11),2)./rho; %cbv in dl^blood/g^tissue
	dhb_content=cbv.*hb.*E0.*100; %deoxyhaemoglobin content in g^dHb/100 g^tissue
	scatter(const.*sum(Vf(:,6:11),2).*(1-Yv),paramsASEtot(1,:),[],dhb_content,'filled');
	hold on;
	plot([0 30],[0 30],'k');
	axis square;
	box on;
	grid on;
	colorbar;
	colormap summer;
	title('Fig. S6a. Apparent R_2^\prime vs SDR predicted R_2^\prime')
	ylabel('Apparent R_2^\prime')
	xlabel('SDR predicted R_2^\prime')


	%FIGURE S6B
	figure;
	scatter(sum(Vf(:,6:11),2),paramsASEtot(2,:),[],E0,'filled');
	hold on;
	plot([0 0.08],[0 0.08],'k');
	axis square;
	box on;
	grid on;
	colorbar;
	colormap autumn;
	title('Fig. S6b. Apparent DBV vs true DBV')
	ylabel('Apparent DBV')
	xlabel('True DBV')
	
	%FIGURE S6C
	figure;
	scatter(E0,paramsASEtot(3,:),[],sum(Vf(:,6:11),2),'filled');
	axis square;
	box on;
	grid on;
	colorbar;
	colormap winter;
	title('Fig. S6c. Apparent OEF vs true OEF')
	ylabel('Apparent OEF')
	xlabel('True OEF')	
	
	%FIGURE S7
	figure;
	scatter(E0,paramsASEtot(2,:)./sum(Vf(:,6:11),2)',[],sum(Vf(:,6:11),2),'filled');        
	box on;
	grid on;
	axis square;
	colorbar;
	colormap winter;
	ylabel('Error in apparent DBV (%)')
	xlabel('True OEF') 
	title('Fig. S7. Error in apparent DBV vs true OEF')	
