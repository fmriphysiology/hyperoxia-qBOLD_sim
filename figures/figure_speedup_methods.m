function figure_speedup_methods(simdir)

	lc=lines(6);

	%simulating scaling for oxygenation
	Rs=[5 500];
	
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0/simvessim_res' num2str(Rs(k))]);
		[sigASEY60(:,k) tauASE]=generate_signal(p,spp,'display',false,'seq','ASE','includeIV',false,'T2EV',Inf);
		load([simdir 'single_vessel_radius_D1-0Y80/simvessim_res' num2str(Rs(k))]);
		[sigASEY80(:,k) tauASE]=generate_signal(p,spp,'display',false,'seq','ASE','includeIV',false,'T2EV',Inf,'Y',0.6);		
	end	
	
	figure;
	hold on;
	for k=1:length(Rs)
		plot(tauASE.*1000,sigASEY60(:,k),'o','color',lc(k,:));
		plot(tauASE.*1000,sigASEY80(:,k),'-','color',lc(k,:));
	end
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([0.6 1]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[0.6:0.1:1]);
	%legend('S_{EV}(R_c=1000\mum)','S_{EV}(R_c=50\mum)','S_{EV}(R_c=10\mum)','S_{EV}(R_c=5\mum)','location','south')
	grid on;
	axis square;
	title('Scaling for oxygenation');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');
	
	%simulating scaled shape function
	Rs=[5 500];
	
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0/simvessim_res' num2str(Rs(k))]);
		[sigASEVf5(:,k) tauASE]=generate_signal(p,spp,'display',false,'seq','ASE','includeIV',false,'T2EV',Inf);
		load([simdir 'single_vessel_radius_D1-0Vf1pc/simvessim_res' num2str(Rs(k))]);
		[sigASEVf1(:,k) tauASE]=generate_signal(p,spp,'display',false,'seq','ASE','includeIV',false,'T2EV',Inf);		
	end	
	
	figure;
	hold on;
	for k=1:length(Rs)
		plot(tauASE.*1000,sigASEVf5(:,k),'o','color',lc(k,:));
		plot(tauASE.*1000,exp(0.05.*log(sigASEVf1(:,k))./0.01),'-','color',lc(k,:));
	end
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([0.6 1]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[0.6:0.1:1]);
	%legend('S_{EV}(R_c=1000\mum)','S_{EV}(R_c=50\mum)','S_{EV}(R_c=10\mum)','S_{EV}(R_c=5\mum)','location','south')
	grid on;
	axis square;
	title('Shape functions');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');
	
	%simulating different levels of diffusion
	RD1=[5 500];
	RD2=[7 700];
	
	for k=1:length(RD1)
		load([simdir 'single_vessel_radius_D1-0/simvessim_res' num2str(RD1(k))]);
		[sigASED1(:,k) tauASE]=generate_signal(p,spp,'display',false,'seq','ASE','includeIV',false,'T2EV',Inf);
		load([simdir 'single_vessel_radius_D2-0/simvessim_res' num2str(RD2(k))]);
		[sigASED2(:,k) tauASE]=generate_signal(p,spp,'display',false,'seq','ASE','includeIV',false,'T2EV',Inf);		
	end
	
	figure;
	hold on;
	for k=1:length(RD1)
		plot(tauASE.*1000,sigASED1(:,k),'o','color',lc(k,:));
		plot(tauASE.*1000,sigASED2(:,k),'-','color',lc(k,:));
	end
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([0.6 1]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[0.6:0.1:1]);
	%legend('S_{EV}(R_c=1000\mum)','S_{EV}(R_c=50\mum)','S_{EV}(R_c=10\mum)','S_{EV}(R_c=5\mum)','location','south')
	grid on;
	axis square;
	title('Diffusion');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');


