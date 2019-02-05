function T=figure_compare_gesse_ase(simdir)

	Rs=[1000 50 10 5];
	Y=0.6;
	Vf=0.03;
	
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASE(:,k) tauASE sigASEev(:,k) sigASEiv(:,k)]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','includeIV',true,'T2EV',Inf,'T2b0',Inf);
		[sigGESSE(:,k) tauGESSE sigGESSEev(:,k) sigGESSEiv(:,k)]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','GESSE','includeIV',true,'T2EV',Inf,'T2b0',Inf);	
	end
	
	%FIGURE 2A
	figure;
	hold on;
	for k=1:length(Rs)
		plot(tauASE.*1000,sigASEev(:,k),'-');
	end
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([0.8 1]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[0.8:0.05:1]);
	legend('S_{EV}(R_c=1000\mum)','S_{EV}(R_c=50\mum)','S_{EV}(R_c=10\mum)','S_{EV}(R_c=5\mum)','location','south')
	grid on;
	title('Fig. 2a. Extravascular ASE signal decay');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');	
	
	%FIGURE 2B
	figure;
	hold on;
	plot(tauASE.*1000,sigASEiv(:,k),'-');
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([0 0.6]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[0:0.1:0.6]);
	legend('S_{IV}','location','south')
	grid on;
	title('Fig. 2b. Intravascular ASE signal decay');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');
	%SUPPLEMENTARY FIGURE

	%FIGURE S4A
	figure;
	hold on;
	for k=1:length(Rs)
		plot(tauGESSE.*1000,sigGESSEev(:,k),'-');
	end
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([0.8 1]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[0.8:0.05:1]);
	legend('S_{EV}(R_c=1000\mum)','S_{EV}(R_c=50\mum)','S_{EV}(R_c=10\mum)','S_{EV}(R_c=5\mum)','location','south')
	grid on;
	title('Fig. S4a. Extravascular GESSE signal decay');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');

	%FIGURE S4B
	figure;
	hold on;
	plot(tauGESSE.*1000,sigGESSEiv(:,k),'-');
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([0 0.6]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[0:0.1:0.6]);
	legend('S_{IV}','location','south')
	grid on;
	title('Fig. S4b. Intravascular GESSE signal decay');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');
	