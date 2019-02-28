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
	deltaW=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(1-Y));
	
	%FIGURE 3A
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
	title('Fig. 3a. Effect of blood volume on apparent R_2^\prime');
	xlabel('Vessel radius (\mum)');
	ylabel('Reversible relaxation rate, R_2^\prime (s^{-1})');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;
	
	%FIGURE 3B	
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
	title('Fig. 3b. Effect of blood volume on apparent DBV');
	xlabel('Vessel radius (\mum)');
	ylabel('Deoxygenated blood volume (%)');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;
	
	%FIGURE 3C
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
	title('Fig. 3c. Effect of blood volume on apparent OEF');
	xlabel('Vessel radius (\mum)');
	ylabel('Oxygen extraction fraction (%)');
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
	
	%FIGURE 3D
	%plot effect on R2'
	figure;
	errorbar(Rs,paramsASED1V3Y60(1,:),paramsASED1V3Y60sd(1,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,paramsASED1V3Y40(1,:),paramsASED1V3Y40sd(1,:),'o','color',lc(2,:));
	errorbar(Rs,paramsASED1V3Y80(1,:),paramsASED1V3Y80sd(1,:),'o','color',lc(3,:));
	R2p=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(1-0.6).*Vf);
	plot(Rs,ones(size(Rs)).*R2p,'--','color',lc(1,:));
	R2p=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(1-0.40).*Vf);
	plot(Rs,ones(size(Rs)).*R2p,'--','color',lc(2,:));
	R2p=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(1-0.8).*Vf);
	plot(Rs,ones(size(Rs)).*R2p,'--','color',lc(3,:));	
	set(gca,'xscale','log');
	xlim([1 1000]);
	title('Fig. 3d. Effect of blood oxygenation on apparent R_2^\prime');
	xlabel('Vessel radius (\mum)');
	ylabel('Reversible relaxation rate, R_2^\prime (s^{-1})');
	grid;
	axis square;
	legend('Y=60%','Y=40%','Y=80%','location','northwest');
	drawnow;
	
	%FIGURE 3E
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
	title('Fig. 3e. Effect of blood oxygenation on apparent DBV');
	xlabel('Vessel radius (\mum)');
	ylabel('Deoxygenated blood volume (%)');
	xlim([1 1000]);
	ylim([0 12]);
	grid;
	axis square;
	legend('Y=60%','Y=40%','Y=80%','location','northwest');
	drawnow;
	
	%FIGURE 3F
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
	title('Fig. 3f. Effect of blood oxygenation on apparent OEF');
	xlabel('Vessel radius (\mum)');
	ylabel('Oxygen extraction fraction (%)');
	xlim([1 1000]);
	ylim([0 100]);
	grid;
	axis square;
	legend('Y=60%','Y=40%','Y=80%','location','northwest');
	drawnow;

	%examine scaling of volume

	%FIGURE 4A
	lc=lines(6);
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
	title('Fig. 4a. Overestimation of DBV');
	xlabel('Vessel radius (\mum)');
	ylabel('Overestimation (%)');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;
	
	%FIGURE 4B
	lc=lc(4:6,:);
	figure;
	errorbar(Rs,100.*paramsASED1V3Y60(2,:)./0.03-100,100.*paramsASED1V3Y60sd(2,:)./0.03,'o','color',lc(1,:));
	hold on;
	errorbar(Rs,100.*paramsASED1V3Y40(2,:)./0.03-100,100.*paramsASED1V3Y40sd(2,:)./0.03,'o','color',lc(2,:));
	errorbar(Rs,100.*paramsASED1V3Y80(2,:)./0.03-100,100.*paramsASED1V3Y80sd(2,:)./0.03,'o','color',lc(3,:));
	set(gca,'xscale','log');
	xlim([1 1000]);
	%ylim([0 12]);
	grid;
	axis square;	
	title('Fig. 4b. Overestimation of DBV');
	xlabel('Vessel radius (\mum)');
	ylabel('Overestimation (%)');
	legend('Y=60%','Y=40%','Y=80%','location','northwest');
	drawnow;
	
	%repeat simulations with intravascular signal included
	
	%change DBV, fix OEF
	
	%ASE signal with D=1, Vf=3, Y=60
	Vf=0.03;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V3Y60f(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',false,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V3Y60f(:,k) paramsASED1V3Y60fsd(:,k)]=calc_qbold_params(p,sigASED1V3Y60f(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end
	fprintf('\n')
	
	%ASE signal with D=1, Vf=5, Y=60
	Vf=0.05;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0/simvessim_res' num2str(Rs(k))]);
		[sigASED1V5Y60f(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',false,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V5Y60f(:,k) paramsASED1V5Y60fsd(:,k)]=calc_qbold_params(p,sigASED1V5Y60f(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end	
	fprintf('\n')
	
	%ASE signal with D=1, Vf=1, Y=60
	Vf=0.01;
	Y=0.6;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf1pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V1Y60f(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',false,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V1Y60f(:,k) paramsASED1V1Y60fsd(:,k)]=calc_qbold_params(p,sigASED1V1Y60f(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end
	fprintf('\n')
	
	lc=lines(6);
	
	%FIGURE 5A
	%plot relative error in R2'
	figure;
	plot(Rs,(paramsASED1V3Y60f(1,:)-paramsASED1V3Y60(1,:))./paramsASED1V3Y60(1,:).*100,'o','color',lc(1,:));
	hold on;
	plot(Rs,(paramsASED1V1Y60f(1,:)-paramsASED1V1Y60(1,:))./paramsASED1V1Y60(1,:).*100,'o','color',lc(2,:));
	plot(Rs,(paramsASED1V5Y60f(1,:)-paramsASED1V5Y60(1,:))./paramsASED1V5Y60(1,:).*100,'o','color',lc(3,:));
	set(gca,'xscale','log');
	xlim([1 1000]);
	ylim([-10 4])
	grid;
	axis square;
	title('Fig. 5a. Effect of intravascular signal on estimate of R_2^\prime');
	xlabel('Vessel radius (\mum)');
	%ylabel('Reversible relaxation rate, R_2^\prime (s^{-1})');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northeast');
	drawnow;
	
	%FIGURE 5B
	%plot effect on DBV
	figure;
	plot(Rs,(paramsASED1V3Y60f(2,:)-paramsASED1V3Y60(2,:))./paramsASED1V3Y60(2,:).*100,'o','color',lc(1,:));
	hold on;
	plot(Rs,(paramsASED1V1Y60f(2,:)-paramsASED1V1Y60(2,:))./paramsASED1V1Y60(2,:).*100,'o','color',lc(2,:));
	plot(Rs,(paramsASED1V5Y60f(2,:)-paramsASED1V5Y60(2,:))./paramsASED1V5Y60(2,:).*100,'o','color',lc(3,:));	
	set(gca,'xscale','log');
	xlim([1 1000]);
	ylim([-10 4])
	grid;
	axis square;	
	title('Fig. 5b. Effect of intravascular signal on estimate of DBV');
	xlabel('Vessel radius (\mum)');
	%ylabel('Deoxygenated blood volume (%)');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northeast');
	drawnow;
	
	%FIGURE 5C
	%plot effect on OEF
	figure;
	plot(Rs,(paramsASED1V3Y60f(3,:)-paramsASED1V3Y60(3,:))./paramsASED1V3Y60(3,:).*100,'o','color',lc(1,:));
	hold on;
	plot(Rs,(paramsASED1V1Y60f(3,:)-paramsASED1V1Y60(3,:))./paramsASED1V1Y60(3,:).*100,'o','color',lc(2,:));
	plot(Rs,(paramsASED1V5Y60f(3,:)-paramsASED1V5Y60(3,:))./paramsASED1V5Y60(3,:).*100,'o','color',lc(3,:));
	set(gca,'xscale','log');
	xlim([1 1000]);
	ylim([-10 4])
	grid;
	axis square;	
	title('Fig. 5c. Effect of intravascular signal on estimate of OEF');
	xlabel('Vessel radius (\mum)');
	%ylabel('Oxygen extraction fraction (%)');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northeast');
	drawnow;
	
	%change OEF, fix DBV
	
	%ASE signal with D=1, Vf=3, Y=40
	Vf=0.03;
	Y=0.4;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V3Y40f(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',false,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V3Y40f(:,k) paramsASED1V3Y40fsd(:,k)]=calc_qbold_params(p,sigASED1V3Y40f(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end	
	fprintf('\n')
	
	%ASE signal with D=1, Vf=3, Y=80
	Vf=0.03;
	Y=0.8;
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V3Y80f(:,k) tauASE]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE,'tau',taus,'includeIV',false,'T2EV',80e-3,'T2b0',189e-3);
		[paramsASED1V3Y80f(:,k) paramsASED1V3Y80fsd(:,k)]=calc_qbold_params(p,sigASED1V3Y80f(:,k),tauASE,tau_cutoff);
		fprintf('.')
	end
	fprintf('\n')	

	lc=lines(6);
	lc=lc(4:6,:);

	%FIGURE 5D
	%plot relative error in R2'
	figure;
	plot(Rs,(paramsASED1V3Y60f(1,:)-paramsASED1V3Y60(1,:))./paramsASED1V3Y60(1,:).*100,'o','color',lc(1,:));
	hold on;
	plot(Rs,(paramsASED1V3Y40f(1,:)-paramsASED1V3Y40(1,:))./paramsASED1V3Y40(1,:).*100,'o','color',lc(2,:));
	plot(Rs,(paramsASED1V3Y80f(1,:)-paramsASED1V3Y80(1,:))./paramsASED1V3Y80(1,:).*100,'o','color',lc(3,:));
	set(gca,'xscale','log');
	xlim([1 1000]);
	ylim([-10 4])
	grid;
	axis square;
	title('Fig. 5d. Effect of intravascular signal on estimate of R_2^\prime');
	xlabel('Vessel radius (\mum)');
	%ylabel('Reversible relaxation rate, R_2^\prime (s^{-1})');
	legend('Y=60%','Y=40%','Y=80%','location','northeast');
	drawnow;
	
	%FIGURE 5E
	%plot effect on DBV
	figure;
	plot(Rs,(paramsASED1V3Y60f(2,:)-paramsASED1V3Y60(2,:))./paramsASED1V3Y60(2,:).*100,'o','color',lc(1,:));
	hold on;
	plot(Rs,(paramsASED1V3Y40f(2,:)-paramsASED1V3Y40(2,:))./paramsASED1V3Y40(2,:).*100,'o','color',lc(2,:));
	plot(Rs,(paramsASED1V3Y80f(2,:)-paramsASED1V3Y80(2,:))./paramsASED1V3Y80(2,:).*100,'o','color',lc(3,:));	
	set(gca,'xscale','log');
	xlim([1 1000]);
	ylim([-10 4])
	grid;
	axis square;	
	title('Fig. 5e. Effect of intravascular signal on estimate of DBV');
	xlabel('Vessel radius (\mum)');
	%ylabel('Deoxygenated blood volume (%)');
	legend('Y=60%','Y=40%','Y=80%','location','northeast');
	drawnow;
	
	%FIGURE 5F
	%plot effect on OEF
	figure;
	plot(Rs,(paramsASED1V3Y60f(3,:)-paramsASED1V3Y60(3,:))./paramsASED1V3Y60(3,:).*100,'o','color',lc(1,:));
	hold on;
	plot(Rs,(paramsASED1V3Y40f(3,:)-paramsASED1V3Y40(3,:))./paramsASED1V3Y40(3,:).*100,'o','color',lc(2,:));
	plot(Rs,(paramsASED1V3Y80f(3,:)-paramsASED1V3Y80(3,:))./paramsASED1V3Y80(3,:).*100,'o','color',lc(3,:));
	set(gca,'xscale','log');
	xlim([1 1000]);
	ylim([-10 4])
	grid;
	axis square;	
	title('Fig. 5f. Effect of intravascular signal on estimate of OEF');
	xlabel('Vessel radius (\mum)');
	%ylabel('Oxygen extraction fraction (%)');
	legend('Y=60%','Y=40%','Y=80%','location','northeast');
	drawnow;
	
	%investigate origins of volume overestimation
	
	%FIGURE 6B
	figure;
	h=area(Rs,100.*paramsASED1V3Y60(2,:));
	h(1).FaceColor=[0.9 0.9 0.9];
	h(1).EdgeColor=[1 1 1];
	hold on;
	p1=errorbar(Rs,-100.*(paramsASED1V3Y60(4,:)-paramsASED1V3Y60(4,end)),100.*paramsASED1V3Y60sd(4,:),'o','color',[255 127 42]./255);
	p2=errorbar(Rs,100.*(paramsASED1V3Y60(4,:)-paramsASED1V3Y60(4,end)+paramsASED1V3Y60(2,:)),100.*sqrt(paramsASED1V3Y60sd(4,:).^2+paramsASED1V3Y60sd(2,:).^2),'o','color',[113 200 55]./255);
	plot(Rs,ones(size(Rs)).*Vf.*100,'--','color',[113 200 55]./255);
	plot(Rs,zeros(size(Rs)),'--','color',[255 127 42]./255);		
	set(gca,'xscale','log');
	title('Fig. 6b. OEF=40%');
	xlabel('Vessel radius (\mum)');
	ylabel('Contribution to DBV (%)');
	xlim([1 1000]);
	ylim([-4 10]);
	grid;
	axis square;
	legend([p1 p2],{'-ln S_{meas}^S(0)','ln S_{extrap}^L(0)'},'location','southeast');
	drawnow;

	%FIGURE 6A
	figure;
	h=area(Rs,100.*paramsASED1V3Y40(2,:));
	h(1).FaceColor=[0.9 0.9 0.9];
	h(1).EdgeColor=[1 1 1];
	hold on;
	p1=errorbar(Rs,-100.*(paramsASED1V3Y40(4,:)-paramsASED1V3Y40(4,end)),100.*paramsASED1V3Y40sd(4,:),'o','color',[255 127 42]./255);
	p2=errorbar(Rs,100.*(paramsASED1V3Y40(4,:)-paramsASED1V3Y40(4,end)+paramsASED1V3Y40(2,:)),100.*sqrt(paramsASED1V3Y40sd(4,:).^2+paramsASED1V3Y40sd(2,:).^2),'o','color',[113 200 55]./255);
	plot(Rs,ones(size(Rs)).*Vf.*100,'--','color',[113 200 55]./255);
	plot(Rs,zeros(size(Rs)),'--','color',[255 127 42]./255);		
	set(gca,'xscale','log');	
	title('Fig. 6a. OEF=60%');
	xlabel('Vessel radius (\mum)');
	ylabel('Contribution to DBV (%)');
	xlim([1 1000]);
	ylim([-4 10]);
	grid;
	axis square;
	legend([p1 p2],{'-ln S_{meas}^S(0)','ln S_{extrap}^L(0)'},'location','southeast');
	drawnow;
	
	%FIGURE 6C
	figure;
	h=area(Rs,100.*paramsASED1V3Y80(2,:));
	h(1).FaceColor=[0.9 0.9 0.9];
	h(1).EdgeColor=[1 1 1];
	hold on;
	p1=errorbar(Rs,-100.*(paramsASED1V3Y80(4,:)-paramsASED1V3Y80(4,end)),100.*paramsASED1V3Y80sd(4,:),'o','color',[255 127 42]./255);
	p2=errorbar(Rs,100.*(paramsASED1V3Y80(4,:)-paramsASED1V3Y80(4,end)+paramsASED1V3Y80(2,:)),100.*sqrt(paramsASED1V3Y80sd(4,:).^2+paramsASED1V3Y80sd(2,:).^2),'o','color',[113 200 55]./255);
	plot(Rs,ones(size(Rs)).*Vf.*100,'--','color',[113 200 55]./255);
	plot(Rs,zeros(size(Rs)),'--','color',[255 127 42]./255);		
	set(gca,'xscale','log');	
	title('Fig. 6c. OEF=20%');
	xlabel('Vessel radius (\mum)');
	ylabel('Contribution to DBV (%)');
	xlim([1 1000]);
	ylim([-4 10]);
	grid;
	axis square;
	legend([p1 p2],{'-ln S_{meas}^S(0)','ln S_{extrap}^L(0)'},'location','southeast');
	drawnow;
		
	%SUPPLEMENTARY FIGURE
	
	lc=lines(6);
	Y=0.6;
	deltaW=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(1-Y));
	
	%FIGURE S5A
	%plot effect on R2'
	figure;
	errorbar(Rs,paramsASED1V3Y60f(1,:),paramsASED1V3Y60fsd(1,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,paramsASED1V1Y60f(1,:),paramsASED1V1Y60fsd(1,:),'o','color',lc(2,:));
	errorbar(Rs,paramsASED1V5Y60f(1,:),paramsASED1V5Y60fsd(1,:),'o','color',lc(3,:));
	plot(Rs,ones(size(Rs)).*deltaW.*0.03,'--','color',lc(1,:));
	plot(Rs,ones(size(Rs)).*deltaW.*0.01,'--','color',lc(2,:));
	plot(Rs,ones(size(Rs)).*deltaW.*0.05,'--','color',lc(3,:));	
	set(gca,'xscale','log');
	xlim([1 1000]);
	grid;
	axis square;
	title('Fig. S5a. Effect of blood volume on apparent R_2^\prime');
	xlabel('Vessel radius (\mum)');
	ylabel('Reversible relaxation rate, R_2^\prime (s^{-1})');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;
	
	%FIGURE S5B
	%plot effect on DBV
	figure;
	errorbar(Rs,100.*paramsASED1V3Y60f(2,:),100.*paramsASED1V3Y60fsd(2,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,100.*paramsASED1V1Y60f(2,:),100.*paramsASED1V1Y60fsd(2,:),'o','color',lc(2,:));
	errorbar(Rs,100.*paramsASED1V5Y60f(2,:),100.*paramsASED1V5Y60fsd(2,:),'o','color',lc(3,:));
	plot(Rs,100.*ones(size(Rs)).*0.03,'--','color',lc(1,:));
	plot(Rs,100.*ones(size(Rs)).*0.01,'--','color',lc(2,:));
	plot(Rs,100.*ones(size(Rs)).*0.05,'--','color',lc(3,:));	
	set(gca,'xscale','log');
	xlim([1 1000]);
	ylim([0 12]);
	grid;
	axis square;	
	title('Fig. S5b. Effect of blood volume on apparent DBV');
	xlabel('Vessel radius (\mum)');
	ylabel('Deoxygenated blood volume (%)');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;
	
	%FIGURE S5C
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
	title('Fig. S5c. Effect of blood volume on apparent OEF');
	xlabel('Vessel radius (\mum)');
	ylabel('Oxygen extraction fraction (%)');
	legend('DBV=3%','DBV=1%','DBV=5%','location','northwest');
	drawnow;	
	
	lc=lines(6);
	lc=lc(4:6,:);
	%deltaW=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(1-Y));
	
	%FIGURE S5D
	%plot effect on R2'
	figure;
	errorbar(Rs,paramsASED1V3Y60f(1,:),paramsASED1V3Y60fsd(1,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,paramsASED1V3Y40f(1,:),paramsASED1V3Y40fsd(1,:),'o','color',lc(2,:));
	errorbar(Rs,paramsASED1V3Y80f(1,:),paramsASED1V3Y80fsd(1,:),'o','color',lc(3,:));
	R2p=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(1-0.6).*Vf);
	plot(Rs,ones(size(Rs)).*R2p,'--','color',lc(1,:));
	R2p=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(1-0.40).*Vf);
	plot(Rs,ones(size(Rs)).*R2p,'--','color',lc(2,:));
	R2p=(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct.*(1-0.8).*Vf);
	plot(Rs,ones(size(Rs)).*R2p,'--','color',lc(3,:));	
	set(gca,'xscale','log');
	xlim([1 1000]);
	title('Fig. S5d. Effect of blood oxygenation on apparent R_2^\prime');
	xlabel('Vessel radius (\mum)');
	ylabel('Reversible relaxation rate, R_2^\prime (s^{-1})');
	grid;
	axis square;
	legend('Y=60%','Y=40%','Y=80%','location','northwest');
	drawnow;
	
	%FIGURE S5E
	%plot effect on DBV
	figure;
	errorbar(Rs,100.*paramsASED1V3Y60f(2,:),100.*paramsASED1V3Y60fsd(2,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,100.*paramsASED1V3Y40f(2,:),100.*paramsASED1V3Y40fsd(2,:),'o','color',lc(2,:));
	errorbar(Rs,100.*paramsASED1V3Y80f(2,:),100.*paramsASED1V3Y80fsd(2,:),'o','color',lc(3,:));
	plot(Rs,100.*ones(size(Rs)).*Vf,'--','color',lc(1,:));
	plot(Rs,100.*ones(size(Rs)).*Vf,'--','color',lc(2,:));
	plot(Rs,100.*ones(size(Rs)).*Vf,'--','color',lc(3,:));	
	set(gca,'xscale','log');
	title('Fig. S5e. Effect of blood oxygenation on apparent DBV');
	xlabel('Vessel radius (\mum)');
	ylabel('Deoxygenated blood volume (%)');
	xlim([1 1000]);
	ylim([0 12]);
	grid;
	axis square;
	legend('Y=60%','Y=40%','Y=80%','location','northwest');
	drawnow;
	
	%FIGURE S5F
	%plot effect on OEF
	figure;
	errorbar(Rs,100.*paramsASED1V3Y60f(3,:),100.*paramsASED1V3Y60fsd(3,:),'o','color',lc(1,:));
	hold on;
	errorbar(Rs,100.*paramsASED1V3Y40f(3,:),100.*paramsASED1V3Y40fsd(3,:),'o','color',lc(2,:));
	errorbar(Rs,100.*paramsASED1V3Y80f(3,:),100.*paramsASED1V3Y80fsd(3,:),'o','color',lc(3,:));
	plot(Rs,100.*ones(size(Rs)).*(1-0.6),'--','color',lc(1,:));
	plot(Rs,100.*ones(size(Rs)).*(1-0.4),'--','color',lc(2,:));
	plot(Rs,100.*ones(size(Rs)).*(1-0.8),'--','color',lc(3,:));	
	set(gca,'xscale','log');
	title('Fig. S5f. Effect of blood oxygenation on apparent OEF');
	xlabel('Vessel radius (\mum)');
	ylabel('Oxygen extraction fraction (%)');
	xlim([1 1000]);
	ylim([0 100]);
	grid;
	axis square;
	legend('Y=60%','Y=40%','Y=80%','location','northwest');
	drawnow;	
	