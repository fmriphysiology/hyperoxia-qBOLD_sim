function figure_speedup_methods(simdir)

	lc=lines(6);

	%simulating scaling for oxygenation
	Rs=[5 500];
	
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0/simvessim_res' num2str(Rs(k))]);
		[sigASEY60(:,k) tauASE sigASEevY60(:,k)]=generate_signal(p,spp,'display',false,'seq','ASE','includeIV',false,'T2EV',Inf);
		load([simdir 'single_vessel_radius_D1-0Y80/simvessim_res' num2str(Rs(k))]);
		[sigASEY80(:,k) tauASE sigASEevY80(:,k)]=generate_signal(p,spp,'display',false,'seq','ASE','includeIV',false,'T2EV',Inf,'Y',0.6);		
	end	
	
	%FIGURE S1a
	figure;
	hold on;
	for k=1:length(Rs)
		plot(tauASE.*1000,sigASEevY60(:,k),'o','color',lc(k,:));
		plot(tauASE.*1000,sigASEevY80(:,k),'-','color',lc(k,:));
	end
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([0.6 1]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[0.6:0.1:1]);
	%legend('S_{EV}(R_c=1000\mum)','S_{EV}(R_c=50\mum)','S_{EV}(R_c=10\mum)','S_{EV}(R_c=5\mum)','location','south')
	grid on;
	axis square;
	title('Fig. S1a. Scaling for different oxygenation levels');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');
	
	%FIGURE S1b
	figure;
	hold on;
	for k=1:length(Rs)
		plot(tauASE.*1000,(sigASEevY60(:,k)-sigASEevY80(:,k))./sigASEevY60(:,k).*100,'x','color',lc(k,:));
	end
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([-3 3]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[-3:1:3]);
	grid on;
	axis square;
	title('Fig. S1b. Percentage error when scaling for different oxygenation levels');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Error (%)');
	
	%simulating scaled shape function
	Rs=[5 500];
	
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0/simvessim_res' num2str(Rs(k))]);
		[sigASEVf5(:,k) tauASE sigASEevVf5(:,k)]=generate_signal(p,spp,'display',false,'seq','ASE','includeIV',false,'T2EV',Inf);
		load([simdir 'single_vessel_radius_D1-0Vf1pc/simvessim_res' num2str(Rs(k))]);
		[sigASEVf1(:,k) tauASE sigASEevVf1(:,k)]=generate_signal(p,spp,'display',false,'seq','ASE','includeIV',false,'T2EV',Inf);		
	end	
	
	%FIGURE S2a
	figure;
	hold on;
	for k=1:length(Rs)
		plot(tauASE.*1000,sigASEevVf5(:,k),'o','color',lc(k,:));
		plot(tauASE.*1000,exp(0.05.*log(sigASEevVf1(:,k))./0.01),'-','color',lc(k,:));
	end
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([0.6 1]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[0.6:0.1:1]);
	%legend('S_{EV}(R_c=1000\mum)','S_{EV}(R_c=50\mum)','S_{EV}(R_c=10\mum)','S_{EV}(R_c=5\mum)','location','south')
	grid on;
	axis square;
	title('Fig. S2a. Scaling for different volume fractions');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');
	
	%FIGURE S2b
	figure;
	hold on;
	for k=1:length(Rs)
		plot(tauASE.*1000,(sigASEevVf5(:,k)-exp(0.05.*log(sigASEevVf1(:,k))./0.01))./sigASEevVf5(:,k).*100,'x','color',lc(k,:));
	end
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([-3 3]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[-3:1:3]);
	grid on;
	axis square;
	title('Fig. S2b. Percentage error when scaling for different volume fractions');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Error (%)');
	
	%simulating combining vessel scales	
	load([simdir 'multi_vessel_radius_speedup/simvessim_res3-7.mat']);
	[sigASERC tauASE sigASEevRC]=generate_signal(p,spp,'display',false,'seq','ASE','includeIV',false,'T2EV',Inf);
	load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res3.mat']);
	[sigASERC3 tauASE sigASEevRC3]=generate_signal(p,spp,'display',false,'seq','ASE','Vf',0.025,'includeIV',false,'T2EV',Inf);
	load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res7.mat']);
	[sigASERC7 tauASE sigASEevRC7]=generate_signal(p,spp,'display',false,'seq','ASE','Vf',0.025,'includeIV',false,'T2EV',Inf);
	
	keyboard;
	
	%FIGURE S3a
	figure;
	hold on;
	plot(tauASE.*1000,sigASEevRC,'o','color',lc(1,:));
	plot(tauASE.*1000,prod([sigASEevRC3 sigASEevRC7],2),'-','color',lc(1,:));
	plot(tauASE.*1000,sigASEevRC3,'-','color',lc(2,:));
	plot(tauASE.*1000,sigASEevRC7,'-','color',lc(3,:));
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([0.6 1]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[0.6:0.1:1]);
	%legend('S_{EV}(R_c=1000\mum)','S_{EV}(R_c=50\mum)','S_{EV}(R_c=10\mum)','S_{EV}(R_c=5\mum)','location','south')
	grid on;
	axis square;
	title('Fig. S3a. Combining multiple vessel radii');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');

	%FIGURE S3b
	figure;
	hold on;
	plot(tauASE.*1000,(sigASEevRC-prod([sigASEevRC3 sigASEevRC7],2))./sigASEevRC.*100,'x','color',lc(1,:));
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([-3 3]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[-3:1:3]);
	%legend('S_{EV}(R_c=1000\mum)','S_{EV}(R_c=50\mum)','S_{EV}(R_c=10\mum)','S_{EV}(R_c=5\mum)','location','south')
	grid on;
	axis square;
	title('Fig. S3b. Percentage error when combining multiple vessel radii');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Error (%)');
	
	%simulating different levels of diffusion
	RD1=[5 500];
	RD2=[7 700];
	
	for k=1:length(RD1)
		load([simdir 'single_vessel_radius_D1-0/simvessim_res' num2str(RD1(k))]);
		[sigASED1(:,k) tauASE sigASEevD1(:,k)]=generate_signal(p,spp,'display',false,'seq','ASE','includeIV',false,'T2EV',Inf);
		load([simdir 'single_vessel_radius_D2-0/simvessim_res' num2str(RD2(k))]);
		[sigASED2(:,k) tauASE sigASEevD2(:,k)]=generate_signal(p,spp,'display',false,'seq','ASE','includeIV',false,'T2EV',Inf);		
	end
	
	%FIGURE S4a 
	figure;
	hold on;
	for k=1:length(RD1)
		plot(tauASE.*1000,sigASEevD1(:,k),'o','color',lc(k,:));
		plot(tauASE.*1000,sigASEevD2(:,k),'-','color',lc(k,:));
	end
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([0.6 1]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[0.6:0.1:1]);
	%legend('S_{EV}(R_c=1000\mum)','S_{EV}(R_c=50\mum)','S_{EV}(R_c=10\mum)','S_{EV}(R_c=5\mum)','location','south')
	grid on;
	axis square;
	title('Fig. S4a. Scaling for different rates of diffusion');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');

	%FIGURE S4b 
	figure;
	hold on;
	for k=1:length(RD1)
		plot(tauASE.*1000,(sigASEevD1(:,k)-sigASEevD2(:,k))./sigASEevD1(:,k).*100,'x','color',lc(k,:));
	end
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([-3 3]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[-3:1:3]);
	grid on;
	axis square;
	title('Fig. S4b. Percentage error when scaling for different rates of diffusion');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Error (%)');
	