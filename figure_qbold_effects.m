function figure_qbold_effects(simdir)

	Rs=[1:10 20:10:100 200:100:1000];
	TE=80e-3;
	tau_cutoff=15e-3; %30e-3; %default - 15e-3
	
	%change DBV, fix OEF
	
	%ASE signal with D=1, Vf=3, Y=60
	Vf=0.03;
	Y=0.6;
	for k=1:length(Rs)
		disp(['Step ' num2str(k) ' of ' num2str(length(Rs))]);
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V3Y60(:,:,k) tauASE paramsASED1V3Y60(:,:,k)]=qbold_bootstrp(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',[0 16:4:64]./1000,'tau_cutoff',tau_cutoff);
	end
	
	%ASE signal with D=1, Vf=5, Y=60
	Vf=0.05;
	Y=0.6;
	for k=1:length(Rs)
		disp(['Step ' num2str(k) ' of ' num2str(length(Rs))]);
		load([simdir 'single_vessel_radius_D1-0/simvessim_res' num2str(Rs(k))]);
		[sigASED1V5Y60(:,:,k) tauASE paramsASED1V5Y60(:,:,k)]=qbold_bootstrp(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',[0 16:4:64]./1000,'tau_cutoff',tau_cutoff);
	end	
	
	%ASE signal with D=1, Vf=1, Y=60
	Vf=0.01;
	Y=0.6;
	for k=1:length(Rs)
		disp(['Step ' num2str(k) ' of ' num2str(length(Rs))]);
		load([simdir 'single_vessel_radius_D1-0Vf1pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V1Y60(:,:,k) tauASE paramsASED1V1Y60(:,:,k)]=qbold_bootstrp(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',[0 16:4:64]./1000,'tau_cutoff',tau_cutoff);
	end

	lc=lines(6);
	Yoff=0.95;
	Y=0.6
	deltaW=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(Yoff-Y));
	
	%plot effect on R2'
	figure;
	errorbar(Rs,squeeze(mean(paramsASED1V3Y60(2,:,:),2)),squeeze(std(paramsASED1V3Y60(2,:,:),[],2)),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,squeeze(mean(paramsASED1V1Y60(2,:,:),2)),squeeze(std(paramsASED1V1Y60(2,:,:),[],2)),'o','color',lc(2,:));
	errorbar(Rs,squeeze(mean(paramsASED1V5Y60(2,:,:),2)),squeeze(std(paramsASED1V5Y60(2,:,:),[],2)),'o','color',lc(3,:));
	plot(Rs,ones(size(Rs)).*deltaW.*0.03,'--','color',lc(1,:));
	plot(Rs,ones(size(Rs)).*deltaW.*0.01,'--','color',lc(2,:));
	plot(Rs,ones(size(Rs)).*deltaW.*0.05,'--','color',lc(3,:));	
	set(gca,'xscale','log');
	xlim([1 1000]);
	grid;
	axis square;
	title('Effect of blood volume on apparent R_2^\prime');
	xlabel('Vessel radius (\mum)');
	ylabel('Reversible relaxation rate, R_2^\prime (s^{-1})');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;
	
	%plot effect on DBV
	figure;
	errorbar(Rs,squeeze(mean(100.*paramsASED1V3Y60(1,:,:),2)),squeeze(std(100.*paramsASED1V3Y60(1,:,:),[],2)),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,squeeze(mean(100.*paramsASED1V1Y60(1,:,:),2)),squeeze(std(100.*paramsASED1V1Y60(1,:,:),[],2)),'o','color',lc(2,:));
	errorbar(Rs,squeeze(mean(100.*paramsASED1V5Y60(1,:,:),2)),squeeze(std(100.*paramsASED1V5Y60(1,:,:),[],2)),'o','color',lc(3,:));
	plot(Rs,100.*ones(size(Rs)).*0.03,'--','color',lc(1,:));
	plot(Rs,100.*ones(size(Rs)).*0.01,'--','color',lc(2,:));
	plot(Rs,100.*ones(size(Rs)).*0.05,'--','color',lc(3,:));	
	set(gca,'xscale','log');
	xlim([1 1000]);
	ylim([0 9]);
	grid;
	axis square;	
	title('Effect of blood volume on apparent DBV');
	xlabel('Vessel radius (\mum)');
	ylabel('Deoxygenated blood volume (dimensionless)');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;
	
	%plot effect on OEF
	figure;
	errorbar(Rs,squeeze(mean(100.*paramsASED1V3Y60(3,:,:),2)),squeeze(std(100.*paramsASED1V3Y60(3,:,:),[],2)),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,squeeze(mean(100.*paramsASED1V1Y60(3,:,:),2)),squeeze(std(100.*paramsASED1V1Y60(3,:,:),[],2)),'o','color',lc(2,:));
	errorbar(Rs,squeeze(mean(100.*paramsASED1V5Y60(3,:,:),2)),squeeze(std(100.*paramsASED1V5Y60(3,:,:),[],2)),'o','color',lc(3,:));
	plot(Rs,100.*ones(size(Rs)).*(1-Y),'--','color',lc(1,:));
	plot(Rs,100.*ones(size(Rs)).*(1-Y),'--','color',lc(2,:));
	plot(Rs,100.*ones(size(Rs)).*(1-Y),'--','color',lc(3,:));	
	set(gca,'xscale','log');
	xlim([1 1000]);
	ylim([0 100]);
	grid;
	axis square;	
	title('Effect of blood volume on apparent OEF');
	xlabel('Vessel radius (\mum)');
	ylabel('Oxygen extraction fraction (dimensionless)');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;

	%change OEF, fix DBV
	
	%ASE signal with D=1, Vf=3, Y=40
	Vf=0.03;
	Y=0.4;
	for k=1:length(Rs)
		disp(['Step ' num2str(k) ' of ' num2str(length(Rs))]);
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V3Y40(:,:,k) tauASE paramsASED1V3Y40(:,:,k)]=qbold_bootstrp(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',[0 16:4:64]./1000,'tau_cutoff',tau_cutoff);
	end	
	
	%ASE signal with D=1, Vf=3, Y=80
	Vf=0.03;
	Y=0.2;
	for k=1:length(Rs)
		disp(['Step ' num2str(k) ' of ' num2str(length(Rs))]);
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V3Y80(:,:,k) tauASE paramsASED1V3Y80(:,:,k)]=qbold_bootstrp(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',[0 16:4:64]./1000,'tau_cutoff',tau_cutoff);
	end

	lc=lines(6);
	lc=lc(4:6,:);
	Yoff=0.95;
	%deltaW=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(Yoff-Y));
	
	%plot effect on R2'
	figure;
	errorbar(Rs,squeeze(mean(paramsASED1V3Y60(2,:,:),2)),squeeze(std(paramsASED1V3Y60(2,:,:),[],2)),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,squeeze(mean(paramsASED1V3Y40(2,:,:),2)),squeeze(std(paramsASED1V3Y40(2,:,:),[],2)),'o','color',lc(2,:));
	errorbar(Rs,squeeze(mean(paramsASED1V3Y80(2,:,:),2)),squeeze(std(paramsASED1V3Y80(2,:,:),[],2)),'o','color',lc(3,:));
	R2p=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(Yoff-0.6).*Vf);
	plot(Rs,ones(size(Rs)).*R2p,'--','color',lc(1,:));
	R2p=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(Yoff-0.40).*Vf);
	plot(Rs,ones(size(Rs)).*R2p,'--','color',lc(2,:));
	R2p=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(Yoff-0.8).*Vf);
	plot(Rs,ones(size(Rs)).*R2p,'--','color',lc(3,:));	
	set(gca,'xscale','log');
	xlim([1 1000]);
	title('Effect of blood oxygenation on apparent R_2^\prime');
	xlabel('Vessel radius (\mum)');
	ylabel('Reversible relaxation rate, R_2^\prime (s^{-1})');
	grid;
	axis square;
	legend('Y=60%','Y=40%','Y=80%','location','northwest');
	drawnow;
	
	%plot effect on DBV
	figure;
	errorbar(Rs,squeeze(mean(100.*paramsASED1V3Y60(1,:,:),2)),squeeze(std(100.*paramsASED1V3Y60(1,:,:),[],2)),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,squeeze(mean(100.*paramsASED1V3Y40(1,:,:),2)),squeeze(std(100.*paramsASED1V3Y40(1,:,:),[],2)),'o','color',lc(2,:));
	errorbar(Rs,squeeze(mean(100.*paramsASED1V3Y80(1,:,:),2)),squeeze(std(100.*paramsASED1V3Y80(1,:,:),[],2)),'o','color',lc(3,:));
	plot(Rs,100.*ones(size(Rs)).*Vf,'--','color',lc(1,:));
	plot(Rs,100.*ones(size(Rs)).*Vf,'--','color',lc(2,:));
	plot(Rs,100.*ones(size(Rs)).*Vf,'--','color',lc(3,:));	
	set(gca,'xscale','log');
	title('Effect of blood oxygenation on apparent DBV');
	xlabel('Vessel radius (\mum)');
	ylabel('Deoxygenated blood volume (dimensionless)');
	xlim([1 1000]);
	ylim([0 9]);
	grid;
	axis square;
	legend('Y=60%','Y=40%','Y=80%','location','northwest');
	drawnow;

	%plot effect on OEF
	figure;
	errorbar(Rs,squeeze(mean(100.*paramsASED1V3Y60(3,:,:),2)),squeeze(std(100.*paramsASED1V3Y60(3,:,:),[],2)),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,squeeze(mean(100.*paramsASED1V3Y40(3,:,:),2)),squeeze(std(100.*paramsASED1V3Y40(3,:,:),[],2)),'o','color',lc(2,:));
	errorbar(Rs,squeeze(mean(100.*paramsASED1V3Y80(3,:,:),2)),squeeze(std(100.*paramsASED1V3Y80(3,:,:),[],2)),'o','color',lc(3,:));
	plot(Rs,100.*ones(size(Rs)).*(1-0.6),'--','color',lc(1,:));
	plot(Rs,100.*ones(size(Rs)).*(1-0.4),'--','color',lc(2,:));
	plot(Rs,100.*ones(size(Rs)).*(1-0.8),'--','color',lc(3,:));	
	set(gca,'xscale','log');
	title('Effect of blood oxygenation on apparent OEF');
	xlabel('Vessel radius (\mum)');
	ylabel('Oxygen extraction fraction (dimensionless)');
	xlim([1 1000]);
	ylim([0 100]);
	grid;
	axis square;
	legend('Y=60%','Y=40%','Y=80%','location','northwest');
	drawnow;

	%examine scaling of volume as a function of oxygenation
	
	X=[ones(size(tauASE(tauASE>=tau_cutoff))) -tauASE(tauASE>tau_cutoff)];
	for k=1:length(Rs)
		aY60(:,:,k)=X\log(sigASED1V3Y60(find(tauASE>tau_cutoff),:,k));
		aY40(:,:,k)=X\log(sigASED1V3Y40(find(tauASE>tau_cutoff),:,k));
		aY80(:,:,k)=X\log(sigASED1V3Y80(find(tauASE>tau_cutoff),:,k));
	end
	
	figure;
	errorbar(Rs,squeeze(mean(100.*aY60(1,:,:),2)),squeeze(std(100.*aY60(1,:,:),[],2)),'o','color',[255 127 42]./255);
	hold on;
	errorbar(Rs,squeeze(mean(100.*log(sigASED1V3Y60(1,:,:)),2)),squeeze(std(100.*log(sigASED1V3Y60(1,:,:)),[],2)),'o','color',[113 200 55]./255);
	plot(Rs,squeeze(mean(100.*aY60(1,:,:),2)),'k');
	plot(Rs,squeeze(mean(100.*log(sigASED1V3Y60(1,:,:)),2)),'k');
	set(gca,'xscale','log');
	title('OEF=60%');
	xlabel('Vessel radius (\mum)');
	ylabel('Contribution to DBV (%)');
	xlim([1 1000]);
	ylim([-8 8]);
	grid;
	axis square;
	legend('ln S_{extrap}^L(0)','ln S_{meas}^S(0)','location','northeast');
	drawnow;
	
	figure;
	errorbar(Rs,squeeze(mean(100.*aY40(1,:,:),2)),squeeze(std(100.*aY40(1,:,:),[],2)),'o','color',[255 127 42]./255);
	hold on;
	errorbar(Rs,squeeze(mean(100.*log(sigASED1V3Y40(1,:,:)),2)),squeeze(std(100.*log(sigASED1V3Y40(1,:,:)),[],2)),'o','color',[113 200 55]./255);
	plot(Rs,squeeze(mean(100.*aY40(1,:,:),2)),'k');
	plot(Rs,squeeze(mean(100.*log(sigASED1V3Y40(1,:,:)),2)),'k');
	set(gca,'xscale','log');	
	title('OEF=40%');
	xlabel('Vessel radius (\mum)');
	ylabel('Contribution to DBV (%)');
	xlim([1 1000]);
	ylim([-8 8]);
	grid;
	axis square;
	legend('ln S_{extrap}^L(0)','ln S_{meas}^S(0)','location','northeast');
	drawnow;
	
	figure;
	errorbar(Rs,squeeze(mean(100.*aY80(1,:,:),2)),squeeze(std(100.*aY80(1,:,:),[],2)),'o','color',[255 127 42]./255);
	hold on;
	errorbar(Rs,squeeze(mean(100.*log(sigASED1V3Y80(1,:,:)),2)),squeeze(std(100.*log(sigASED1V3Y80(1,:,:)),[],2)),'o','color',[113 200 55]./255);
	plot(Rs,squeeze(mean(100.*aY80(1,:,:),2)),'k');
	plot(Rs,squeeze(mean(100.*log(sigASED1V3Y80(1,:,:)),2)),'k');
	set(gca,'xscale','log');	
	title('OEF=20%');
	xlabel('Vessel radius (\mum)');
	ylabel('Contribution to DBV (%)');
	xlim([1 1000]);
	ylim([-8 8]);
	grid;
	axis square;
	legend('ln S_{extrap}^L(0)','ln S_{meas}^S(0)','location','northeast');
	drawnow;
	
	keyboard;
	