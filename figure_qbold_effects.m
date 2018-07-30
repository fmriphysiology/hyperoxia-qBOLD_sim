function figure_qbold_effects(simdir)

	Rs=[1:10 20:10:100 200:100:1000];
	TE=80e-3;
	tau_cutoff=15e-3; 
	taus=[0 16:4:64]./1000;
	
	%change DBV, fix OEF
	
	%ASE signal with D=1, Vf=3, Y=60
	Vf=0.03;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V3Y60(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',true,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V3Y60(:,k) paramsASED1V3Y60sd(:,k)]=calc_qbold_params(p,sigASED1V3Y60(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end
	fprintf('\n')
	
	%ASE signal with D=1, Vf=5, Y=60
	Vf=0.05;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0/simvessim_res' num2str(Rs(k))]);
		[sigASED1V5Y60(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',true,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V5Y60(:,k) paramsASED1V5Y60sd(:,k)]=calc_qbold_params(p,sigASED1V5Y60(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end	
	fprintf('\n')
	
	%ASE signal with D=1, Vf=1, Y=60
	Vf=0.01;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf1pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V1Y60(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',true,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V1Y60(:,k) paramsASED1V1Y60sd(:,k)]=calc_qbold_params(p,sigASED1V1Y60(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end
	fprintf('\n')

	lc=lines(6);
	Y=0.6;
	deltaW=(4/3*pi*p.gamma.*7.*p.deltaChi0.*p.Hct.*(1-Y));
	
	%plot effect on R2'
	figure;
	errorbar(Rs,paramsASED1V3Y60(1,:),paramsASED1V3Y60sd(1,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,paramsASED1V1Y60(1,:),paramsASED1V1Y60sd(1,:),'o','color',lc(2,:));
	errorbar(Rs,paramsASED1V5Y60(1,:),paramsASED1V5Y60sd(1,:),'o','color',lc(3,:));
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
	errorbar(Rs,100.*paramsASED1V3Y60(2,:),100.*paramsASED1V3Y60sd(2,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,100.*paramsASED1V1Y60(2,:),100.*paramsASED1V1Y60sd(2,:),'o','color',lc(2,:));
	errorbar(Rs,100.*paramsASED1V5Y60(2,:),100.*paramsASED1V5Y60sd(2,:),'o','color',lc(3,:));
	plot(Rs,100.*ones(size(Rs)).*0.03,'--','color',lc(1,:));
	plot(Rs,100.*ones(size(Rs)).*0.01,'--','color',lc(2,:));
	plot(Rs,100.*ones(size(Rs)).*0.05,'--','color',lc(3,:));	
	set(gca,'xscale','log');
	xlim([1 1000]);
	ylim([0 12]);
	grid;
	axis square;	
	title('Effect of blood volume on apparent DBV');
	xlabel('Vessel radius (\mum)');
	ylabel('Deoxygenated blood volume (%)');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;
	
	%plot effect on OEF
	figure;
	errorbar(Rs,100.*paramsASED1V3Y60(3,:),100.*paramsASED1V3Y60sd(3,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,100.*paramsASED1V1Y60(3,:),100.*paramsASED1V1Y60sd(3,:),'o','color',lc(2,:));
	errorbar(Rs,100.*paramsASED1V5Y60(3,:),100.*paramsASED1V5Y60sd(3,:),'o','color',lc(3,:));
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
	ylabel('Oxygen extraction fraction (%)');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;

	%examine scaling of volume

	figure;
	errorbar(Rs,100.*paramsASED1V3Y60(2,:)./0.03-100,100.*paramsASED1V3Y60sd(2,:)./0.03,'o','color',lc(1,:));
	hold on;
	errorbar(Rs,100.*paramsASED1V1Y60(2,:)./0.01-100,100.*paramsASED1V1Y60sd(2,:)./0.01,'o','color',lc(2,:));
	errorbar(Rs,100.*paramsASED1V5Y60(2,:)./0.05-100,100.*paramsASED1V5Y60sd(2,:)./0.05,'o','color',lc(3,:));
	set(gca,'xscale','log');
	xlim([1 1000]);
	%ylim([0 12]);
	grid;
	axis square;	
	title('Overestimation of DBV');
	xlabel('Vessel radius (\mum)');
	ylabel('Overestimation (%)');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;	

	%change OEF, fix DBV
	
	%ASE signal with D=1, Vf=3, Y=40
	Vf=0.03;
	Y=0.4;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V3Y40(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',true,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V3Y40(:,k) paramsASED1V3Y40sd(:,k)]=calc_qbold_params(p,sigASED1V3Y40(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end	
	fprintf('\n')
	
	%ASE signal with D=1, Vf=3, Y=80
	Vf=0.03;
	Y=0.8;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V3Y80(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',true,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V3Y80(:,k) paramsASED1V3Y80sd(:,k)]=calc_qbold_params(p,sigASED1V3Y80(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end
	fprintf('\n')

	lc=lines(6);
	lc=lc(4:6,:);
	%deltaW=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(1-Y));
	
	%plot effect on R2'
	figure;
	errorbar(Rs,paramsASED1V3Y60(1,:),paramsASED1V3Y60sd(1,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,paramsASED1V3Y40(1,:),paramsASED1V3Y40sd(1,:),'o','color',lc(2,:));
	errorbar(Rs,paramsASED1V3Y80(1,:),paramsASED1V3Y80sd(1,:),'o','color',lc(3,:));
	R2p=(4/3*pi*p.gamma.*7.*p.deltaChi0.*p.Hct.*(1-0.6).*Vf);
	plot(Rs,ones(size(Rs)).*R2p,'--','color',lc(1,:));
	R2p=(4/3*pi*p.gamma.*7.*p.deltaChi0.*p.Hct.*(1-0.40).*Vf);
	plot(Rs,ones(size(Rs)).*R2p,'--','color',lc(2,:));
	R2p=(4/3*pi*p.gamma.*7.*p.deltaChi0.*p.Hct.*(1-0.8).*Vf);
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
	errorbar(Rs,100.*paramsASED1V3Y60(2,:),100.*paramsASED1V3Y60sd(2,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,100.*paramsASED1V3Y40(2,:),100.*paramsASED1V3Y40sd(2,:),'o','color',lc(2,:));
	errorbar(Rs,100.*paramsASED1V3Y80(2,:),100.*paramsASED1V3Y80sd(2,:),'o','color',lc(3,:));
	plot(Rs,100.*ones(size(Rs)).*Vf,'--','color',lc(1,:));
	plot(Rs,100.*ones(size(Rs)).*Vf,'--','color',lc(2,:));
	plot(Rs,100.*ones(size(Rs)).*Vf,'--','color',lc(3,:));	
	set(gca,'xscale','log');
	title('Effect of blood oxygenation on apparent DBV');
	xlabel('Vessel radius (\mum)');
	ylabel('Deoxygenated blood volume (%)');
	xlim([1 1000]);
	ylim([0 12]);
	grid;
	axis square;
	legend('Y=60%','Y=40%','Y=80%','location','northwest');
	drawnow;

	%plot effect on OEF
	figure;
	errorbar(Rs,100.*paramsASED1V3Y60(3,:),100.*paramsASED1V3Y60sd(3,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,100.*paramsASED1V3Y40(3,:),100.*paramsASED1V3Y40sd(3,:),'o','color',lc(2,:));
	errorbar(Rs,100.*paramsASED1V3Y80(3,:),100.*paramsASED1V3Y80sd(3,:),'o','color',lc(3,:));
	plot(Rs,100.*ones(size(Rs)).*(1-0.6),'--','color',lc(1,:));
	plot(Rs,100.*ones(size(Rs)).*(1-0.4),'--','color',lc(2,:));
	plot(Rs,100.*ones(size(Rs)).*(1-0.8),'--','color',lc(3,:));	
	set(gca,'xscale','log');
	title('Effect of blood oxygenation on apparent OEF');
	xlabel('Vessel radius (\mum)');
	ylabel('Oxygen extraction fraction (%)');
	xlim([1 1000]);
	ylim([0 100]);
	grid;
	axis square;
	legend('Y=60%','Y=40%','Y=80%','location','northwest');
	drawnow;

	%examine scaling of volume as a function of oxygenation
	
	figure;
	errorbar(Rs,100.*(paramsASED1V3Y60(4,:)-paramsASED1V3Y60(4,end)),100.*paramsASED1V3Y60sd(4,:),'o','color',[255 127 42]./255);
	hold on;
	errorbar(Rs,100.*(paramsASED1V3Y60(4,:)-paramsASED1V3Y60(4,end)+paramsASED1V3Y60(2,:)),100.*sqrt(paramsASED1V3Y60sd(4,:).^2+paramsASED1V3Y60sd(2,:).^2),'o','color',[113 200 55]./255);
	plot(Rs,100.*(paramsASED1V3Y60(4,:)-paramsASED1V3Y60(4,end)),'k');
	plot(Rs,100.*(paramsASED1V3Y60(4,:)-paramsASED1V3Y60(4,end)+paramsASED1V3Y60(2,:)),'k');
	plot(Rs,ones(size(Rs)).*Vf.*100,'--','color',[113 200 55]./255);
	plot(Rs,zeros(size(Rs)),'--','color',[255 127 42]./255);		
	set(gca,'xscale','log');
	title('OEF=40%');
	xlabel('Vessel radius (\mum)');
	ylabel('Contribution to DBV (%)');
	xlim([1 1000]);
	ylim([-8 6]);
	grid;
	axis square;
	legend('ln S_{meas}^S(0)','ln S_{extrap}^L(0)','location','southeast');
	drawnow;

	figure;
	errorbar(Rs,100.*(paramsASED1V3Y40(4,:)-paramsASED1V3Y40(4,end)),100.*paramsASED1V3Y40sd(4,:),'o','color',[255 127 42]./255);
	hold on;
	errorbar(Rs,100.*(paramsASED1V3Y40(4,:)-paramsASED1V3Y40(4,end)+paramsASED1V3Y40(2,:)),100.*sqrt(paramsASED1V3Y40sd(4,:).^2+paramsASED1V3Y40sd(2,:).^2),'o','color',[113 200 55]./255);
	plot(Rs,100.*(paramsASED1V3Y40(4,:)-paramsASED1V3Y40(4,end)),'k');
	plot(Rs,100.*(paramsASED1V3Y40(4,:)-paramsASED1V3Y40(4,end)+paramsASED1V3Y40(2,:)),'k');
	plot(Rs,ones(size(Rs)).*Vf.*100,'--','color',[113 200 55]./255);
	plot(Rs,zeros(size(Rs)),'--','color',[255 127 42]./255);		
	set(gca,'xscale','log');	
	title('OEF=60%');
	xlabel('Vessel radius (\mum)');
	ylabel('Contribution to DBV (%)');
	xlim([1 1000]);
	ylim([-8 6]);
	grid;
	axis square;
	legend('ln S_{meas}^S(0)','ln S_{extrap}^L(0)','location','southeast');
	drawnow;
	
	figure;
	errorbar(Rs,100.*(paramsASED1V3Y80(4,:)-paramsASED1V3Y80(4,end)),100.*paramsASED1V3Y80sd(4,:),'o','color',[255 127 42]./255);
	hold on;
	errorbar(Rs,100.*(paramsASED1V3Y80(4,:)-paramsASED1V3Y80(4,end)+paramsASED1V3Y80(2,:)),100.*sqrt(paramsASED1V3Y80sd(4,:).^2+paramsASED1V3Y80sd(2,:).^2),'o','color',[113 200 55]./255);
	plot(Rs,100.*(paramsASED1V3Y80(4,:)-paramsASED1V3Y80(4,end)),'k');
	plot(Rs,100.*(paramsASED1V3Y80(4,:)-paramsASED1V3Y80(4,end)+paramsASED1V3Y80(2,:)),'k');
	plot(Rs,ones(size(Rs)).*Vf.*100,'--','color',[113 200 55]./255);
	plot(Rs,zeros(size(Rs)),'--','color',[255 127 42]./255);		
	set(gca,'xscale','log');	
	title('OEF=20%');
	xlabel('Vessel radius (\mum)');
	ylabel('Contribution to DBV (%)');
	xlim([1 1000]);
	ylim([-8 6]);
	grid;
	axis square;
	legend('ln S_{meas}^S(0)','ln S_{extrap}^L(0)','location','southeast');
	drawnow;
	
	%keyboard;
	