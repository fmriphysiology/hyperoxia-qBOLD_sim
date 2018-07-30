function figure_investigate_te(simdir)

	Rs=[1:10 20:10:100 200:100:1000];
	tau_cutoff=15e-3; 
	taus=[0 16:4:32]./1000;
	
	%TE=80ms
	TE=80e-3;
	
	%ASE signal with D=1, Vf=3, Y=60
	Vf=0.03;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V3Y60TE80(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',true,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V3Y60TE80(:,k) paramsASED1V3Y60TE80sd(:,k)]=calc_qbold_params(p,sigASED1V3Y60TE80(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end
	fprintf('\n')
	
	%ASE signal with D=1, Vf=5, Y=60
	Vf=0.05;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0/simvessim_res' num2str(Rs(k))]);
		[sigASED1V5Y60TE80(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',true,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V5Y60TE80(:,k) paramsASED1V5Y60TE80sd(:,k)]=calc_qbold_params(p,sigASED1V5Y60TE80(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end	
	fprintf('\n')
	
	%ASE signal with D=1, Vf=1, Y=60
	Vf=0.01;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf1pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V1Y60TE80(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',true,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V1Y60TE80(:,k) paramsASED1V1Y60TE80sd(:,k)]=calc_qbold_params(p,sigASED1V1Y60TE80(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end
	fprintf('\n')

	lc=lines(6);
	Y=0.6;
	deltaW=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(1-Y));
		
	%plot effect on DBV
	figure;
	errorbar(Rs,100.*paramsASED1V3Y60TE80(2,:),100.*paramsASED1V3Y60TE80sd(2,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,100.*paramsASED1V1Y60TE80(2,:),100.*paramsASED1V1Y60TE80sd(2,:),'o','color',lc(2,:));
	errorbar(Rs,100.*paramsASED1V5Y60TE80(2,:),100.*paramsASED1V5Y60TE80sd(2,:),'o','color',lc(3,:));
	plot(Rs,100.*ones(size(Rs)).*0.03,'--','color',lc(1,:));
	plot(Rs,100.*ones(size(Rs)).*0.01,'--','color',lc(2,:));
	plot(Rs,100.*ones(size(Rs)).*0.05,'--','color',lc(3,:));	
	set(gca,'xscale','log');
	xlim([1 1000]);
	ylim([0 12]);
	grid;
	axis square;	
	title('Effect of echo time on DBV (TE=80ms)');
	xlabel('Vessel radius (\mum)');
	ylabel('Deoxygenated blood volume (%)');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;

	%examine scaling of volume as a function of oxygenation
	Vf=0.03;
	figure;
	errorbar(Rs,100.*(paramsASED1V3Y60TE80(4,:)-paramsASED1V3Y60TE80(4,end)),100.*paramsASED1V3Y60TE80sd(4,:),'o','color',[255 127 42]./255);
	hold on;
	errorbar(Rs,100.*(paramsASED1V3Y60TE80(4,:)-paramsASED1V3Y60TE80(4,end)+paramsASED1V3Y60TE80(2,:)),100.*sqrt(paramsASED1V3Y60TE80sd(4,:).^2+paramsASED1V3Y60TE80sd(2,:).^2),'o','color',[113 200 55]./255);
	plot(Rs,100.*(paramsASED1V3Y60TE80(4,:)-paramsASED1V3Y60TE80(4,end)),'k');
	plot(Rs,100.*(paramsASED1V3Y60TE80(4,:)-paramsASED1V3Y60TE80(4,end)+paramsASED1V3Y60TE80(2,:)),'k');
	plot(Rs,ones(size(Rs)).*Vf.*100,'--','color',[113 200 55]./255);
	plot(Rs,zeros(size(Rs)),'--','color',[255 127 42]./255);		
	set(gca,'xscale','log');
	title('Contributions to DBV (TE=80ms, OEF=40%)');
	xlabel('Vessel radius (\mum)');
	ylabel('Contribution to DBV (%)');
	xlim([1 1000]);
	ylim([-8 6]);
	grid;
	axis square;
	legend('ln S_{meas}^S(0)','ln S_{extrap}^L(0)','location','southeast');
	drawnow;

	%TE=48ms
	TE=48e-3;
	
	%ASE signal with D=1, Vf=3, Y=60
	Vf=0.03;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V3Y60TE48(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',true,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V3Y60TE48(:,k) paramsASED1V3Y60TE48sd(:,k)]=calc_qbold_params(p,sigASED1V3Y60TE48(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end
	fprintf('\n')
	
	%ASE signal with D=1, Vf=5, Y=60
	Vf=0.05;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0/simvessim_res' num2str(Rs(k))]);
		[sigASED1V5Y60TE48(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',true,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V5Y60TE48(:,k) paramsASED1V5Y60TE48sd(:,k)]=calc_qbold_params(p,sigASED1V5Y60TE48(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end	
	fprintf('\n')
	
	%ASE signal with D=1, Vf=1, Y=60
	Vf=0.01;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf1pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V1Y60TE48(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',true,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V1Y60TE48(:,k) paramsASED1V1Y60TE48sd(:,k)]=calc_qbold_params(p,sigASED1V1Y60TE48(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end
	fprintf('\n')

	lc=lines(6);
	Y=0.6;
	deltaW=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(1-Y));
		
	%plot effect on DBV
	figure;
	errorbar(Rs,100.*paramsASED1V3Y60TE48(2,:),100.*paramsASED1V3Y60TE48sd(2,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,100.*paramsASED1V1Y60TE48(2,:),100.*paramsASED1V1Y60TE48sd(2,:),'o','color',lc(2,:));
	errorbar(Rs,100.*paramsASED1V5Y60TE48(2,:),100.*paramsASED1V5Y60TE48sd(2,:),'o','color',lc(3,:));
	plot(Rs,100.*ones(size(Rs)).*0.03,'--','color',lc(1,:));
	plot(Rs,100.*ones(size(Rs)).*0.01,'--','color',lc(2,:));
	plot(Rs,100.*ones(size(Rs)).*0.05,'--','color',lc(3,:));	
	set(gca,'xscale','log');
	xlim([1 1000]);
	ylim([0 12]);
	grid;
	axis square;	
	title('Effect of echo time on DBV (TE=48ms)');
	xlabel('Vessel radius (\mum)');
	ylabel('Deoxygenated blood volume (%)');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;

	%examine scaling of volume as a function of oxygenation
	Vf=0.03;
	figure;
	errorbar(Rs,100.*(paramsASED1V3Y60TE48(4,:)-paramsASED1V3Y60TE48(4,end)),100.*paramsASED1V3Y60TE48sd(4,:),'o','color',[255 127 42]./255);
	hold on;
	errorbar(Rs,100.*(paramsASED1V3Y60TE48(4,:)-paramsASED1V3Y60TE48(4,end)+paramsASED1V3Y60TE48(2,:)),100.*sqrt(paramsASED1V3Y60TE48sd(4,:).^2+paramsASED1V3Y60TE48sd(2,:).^2),'o','color',[113 200 55]./255);
	plot(Rs,100.*(paramsASED1V3Y60TE48(4,:)-paramsASED1V3Y60TE48(4,end)),'k');
	plot(Rs,100.*(paramsASED1V3Y60TE48(4,:)-paramsASED1V3Y60TE48(4,end)+paramsASED1V3Y60TE48(2,:)),'k');
	plot(Rs,ones(size(Rs)).*Vf.*100,'--','color',[113 200 55]./255);
	plot(Rs,zeros(size(Rs)),'--','color',[255 127 42]./255);		
	set(gca,'xscale','log');
	title('Contributions to DBV (TE=48ms, OEF=40%)');
	xlabel('Vessel radius (\mum)');
	ylabel('Contribution to DBV (%)');
	xlim([1 1000]);
	ylim([-8 6]);
	grid;
	axis square;
	legend('ln S_{meas}^S(0)','ln S_{extrap}^L(0)','location','southeast');
	drawnow;

	%TE=112ms
	TE=112e-3;
	
	%ASE signal with D=1, Vf=3, Y=60
	Vf=0.03;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V3Y60TE112(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',true,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V3Y60TE112(:,k) paramsASED1V3Y60TE112sd(:,k)]=calc_qbold_params(p,sigASED1V3Y60TE112(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end
	fprintf('\n')
	
	%ASE signal with D=1, Vf=5, Y=60
	Vf=0.05;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0/simvessim_res' num2str(Rs(k))]);
		[sigASED1V5Y60TE112(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',true,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V5Y60TE112(:,k) paramsASED1V5Y60TE112sd(:,k)]=calc_qbold_params(p,sigASED1V5Y60TE112(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end	
	fprintf('\n')
	
	%ASE signal with D=1, Vf=1, Y=60
	Vf=0.01;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf1pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V1Y60TE112(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',true,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V1Y60TE112(:,k) paramsASED1V1Y60TE112sd(:,k)]=calc_qbold_params(p,sigASED1V1Y60TE112(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end
	fprintf('\n')

	lc=lines(6);
	Y=0.6;
	deltaW=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(1-Y));
		
	%plot effect on DBV
	figure;
	errorbar(Rs,100.*paramsASED1V3Y60TE112(2,:),100.*paramsASED1V3Y60TE112sd(2,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,100.*paramsASED1V1Y60TE112(2,:),100.*paramsASED1V1Y60TE112sd(2,:),'o','color',lc(2,:));
	errorbar(Rs,100.*paramsASED1V5Y60TE112(2,:),100.*paramsASED1V5Y60TE112sd(2,:),'o','color',lc(3,:));
	plot(Rs,100.*ones(size(Rs)).*0.03,'--','color',lc(1,:));
	plot(Rs,100.*ones(size(Rs)).*0.01,'--','color',lc(2,:));
	plot(Rs,100.*ones(size(Rs)).*0.05,'--','color',lc(3,:));	
	set(gca,'xscale','log');
	xlim([1 1000]);
	ylim([0 12]);
	grid;
	axis square;	
	title('Effect of echo time on DBV (TE=112ms)');
	xlabel('Vessel radius (\mum)');
	ylabel('Deoxygenated blood volume (%)');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;

	%examine scaling of volume as a function of oxygenation
	Vf=0.03;
	figure;
	errorbar(Rs,100.*(paramsASED1V3Y60TE112(4,:)-paramsASED1V3Y60TE112(4,end)),100.*paramsASED1V3Y60TE112sd(4,:),'o','color',[255 127 42]./255);
	hold on;
	errorbar(Rs,100.*(paramsASED1V3Y60TE112(4,:)-paramsASED1V3Y60TE112(4,end)+paramsASED1V3Y60TE112(2,:)),100.*sqrt(paramsASED1V3Y60TE112sd(4,:).^2+paramsASED1V3Y60TE112sd(2,:).^2),'o','color',[113 200 55]./255);
	plot(Rs,100.*(paramsASED1V3Y60TE112(4,:)-paramsASED1V3Y60TE112(4,end)),'k');
	plot(Rs,100.*(paramsASED1V3Y60TE112(4,:)-paramsASED1V3Y60TE112(4,end)+paramsASED1V3Y60TE112(2,:)),'k');
	plot(Rs,ones(size(Rs)).*Vf.*100,'--','color',[113 200 55]./255);
	plot(Rs,zeros(size(Rs)),'--','color',[255 127 42]./255);		
	set(gca,'xscale','log');
	title('Contributions to DBV (TE=112ms, OEF=40%)');
	xlabel('Vessel radius (\mum)');
	ylabel('Contribution to DBV (%)');
	xlim([1 1000]);
	ylim([-8 6]);
	grid;
	axis square;
	legend('ln S_{meas}^S(0)','ln S_{extrap}^L(0)','location','southeast');
	drawnow;

keyboard;