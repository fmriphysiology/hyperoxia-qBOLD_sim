function figure_sdr(simdir,bids_dir)

	TE=80e-3;
	tauASE=[-28:4:64]./1000;
	
	Ya=0.98;
	E0=0.4;
	Yv=Ya.*(1-E0);
	k=0.4;
	Yc=Ya*k+Yv*(1-k);
	Y=Yv;
	
	Vtot=0.05.*0.793;
	Vf=Vtot;
	
	[sigASEsdr sigASElong sigASEshort]=SDRqBOLD(Y,Vf,tauASE);
	se=find(tauASE==0);
	sigASEsdrn=sigASEsdr./mean(sigASEsdr(se-1:se+1));	
	sigASElongn=sigASElong./mean(sigASEsdr(se-1:se+1));	
	sigASEshortn=sigASEshort./mean(sigASEsdr(se-1:se+1));	

	if exist([bids_dir '/derivatives/group_results.mat'])
		load([bids_dir '/derivatives/group_results.mat'],'tcn')
	
		p=gentemplate;
		p.Y=Y;
		p.vesselFraction=Vf;
		params=calc_qbold_params(p,mean(tcn,2),tauASE,15e-3);

		Ya=1;
		E0=params(3);
		Yv=Ya.*(1-E0);
		k=0.4;
		Yc=Ya*k+Yv*(1-k);
		Y=Yv;
	
		Vtot=params(2);
		Vf=Vtot;
		
		[sigASEsdr2 sigASElong2 sigASEshort2]=SDRqBOLD(Y,Vf,tauASE);
		se=find(tauASE==0);
		sigASEsdr2n=sigASEsdr2./mean(sigASEsdr2(se-1:se+1));	
		sigASElong2n=sigASElong2./mean(sigASEsdr2(se-1:se+1));	
		sigASEshort2n=sigASEshort2./mean(sigASEsdr2(se-1:se+1));	
		
	end

	lc=lines(6);
		
	figure;
	if exist([bids_dir '/derivatives/group_results.mat'])
		plot(tauASE.*1000,tcn,'color',[0.5 0.5 0.5])
	end
	hold on;
	plot(tauASE.*1000,sigASEsdrn,'color',lc(2,:),'linewidth',3)
	%plot(tauASE.*1000,sigASEshortn,'color',lc(2,:),'linewidth',3)
	%plot(tauASE.*1000,sigASElongn,'color',lc(3,:),'linewidth',3)
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([0.8 1.02]);
	set(gca,'xtick',[-28:14:56]);
	set(gca,'ytick',[0.8:0.05:1]);
	grid on;
	axis square;
	title('Multiple vessel scale simulations: SDR (E_0=0.4, V=0.039)');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');

	figure;
	if exist([bids_dir '/derivatives/group_results.mat'])
		plot(tauASE.*1000,tcn,'color',[0.5 0.5 0.5])
	end
	hold on;
	plot(tauASE.*1000,sigASEsdr2n,'color',lc(2,:),'linewidth',3)
	%plot(tauASE.*1000,sigASEshort2n,'color',lc(2,:),'linewidth',3)
	%plot(tauASE.*1000,sigASElong2n,'color',lc(3,:),'linewidth',3)
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([0.8 1.02]);
	set(gca,'xtick',[-28:14:56]);
	set(gca,'ytick',[0.8:0.05:1]);
	grid on;
	axis square;
	title('Multiple vessel scale simulations: SDR (E_0=0.2, V=0.048)');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');
		
	function [sigASEsdr sigASElong sigASEshort]=SDRqBOLD(Y,Vf,tauASE)
	
	deltaChi0=0.264e-6;
	Hct=0.4;
	B0=3;
	gyro=2.6754e+08;

	deltaW=(gyro).*4./3.*pi.*deltaChi0.*Hct.*(1-Y).*B0;
	tc=1./deltaW;
	
	%r2p in long tau regime
	r2plc=Vf.*deltaW; %for cylinders

	%long tau signal
	sigASElong=exp(-r2plc.*abs(tauASE)+Vf); %for cylinders +ve

	%short tau signal
	sigASEshort=exp(-0.3.*Vf.*(deltaW.*tauASE).^2); %for cylinders

	for k=1:length(tauASE)

		eqn=@(u)(2+u).*sqrt(1-u).*(1-besselj(0,1.5.*abs(tauASE(k)).*deltaW.*u))./u.^2;
		fc(k,:)=quad(eqn,0,1)./3;

	end

	%complete signal at all tau
	sigASEsdr=exp(-Vf.*fc); %for random cylinders
	
	return;