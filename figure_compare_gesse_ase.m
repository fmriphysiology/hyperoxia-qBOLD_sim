function T=figure_compare_gesse_ase(simdir)

	Rs=[1000 50 10 5];
	Y=0.6;
	Vf=0.03;
	
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASE(:,k) tauASE sigASEev(:,k) sigASEiv(:,k)]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','includeIV',true,'T2EV',Inf,'T2b0',Inf);
		[sigGESSE(:,k) tauGESSE sigGESSEev(:,k) sigGESSEiv(:,k)]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','GESSE','includeIV',true,'T2EV',Inf,'T2b0',Inf);	
	end
	
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
	title('Extravascular GESSE signal decay');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');
	
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
	title('Extravascular ASE signal decay');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');	

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
	title('Intravascular GESSE signal decay');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');
	
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
	title('Intravascular ASE signal decay');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');

	RsW=1000;
	load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(RsW)]);
	[sigASEW tauASE sigASEevW sigASEivW]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','includeIV',true,'T2EV',80e-3,'T2b0',189e-3);

	figure;
	hold on;
	plot(tauASE.*1000,sigASEW,'k-');
	plot(tauASE.*1000,sigASEevW.*(1-Vf),'-');
	plot(tauASE.*1000,sigASEivW.*Vf,'-');
	box;
	xlim([min(tauASE.*1000) max(tauASE.*1000)]);
	ylim([0 0.5]);
	set(gca,'xtick',[-60:30:60]);
	set(gca,'ytick',[0:0.1:0.5]);
	legend('S_{TOT}','S_{EV}(R_c=1000\mum)','S_{IV}','location','east')
	grid on;
	title('Effect of intravascular signal on ASE qBOLD');
	xlabel('Spin echo displacement time, \tau (ms)');
	ylabel('Signal fraction (arb.)');
	
	[paramsASEW paramsASEW_SD]=calc_qbold_params(p,sigASEW,tauASE,15e-3);
	[paramsASEevW paramsASEevW_SD]=calc_qbold_params(p,sigASEevW,tauASE,15e-3);
	
	Variables={'R2prime' 'DBV' 'OEF'}';
	WithIV=paramsASEW(1:3);
	WithIV_SD=paramsASEW_SD(1:3);
	WithoutIV=paramsASEevW(1:3);
	WithoutIV_SD=paramsASEevW_SD(1:3);
	[h(1) Ztest_pval(1,:)]=ztest(WithIV(1)-WithoutIV(1),0,(WithIV_SD(1)+WithoutIV_SD(1))/2);
	[h(2) Ztest_pval(2,:)]=ztest(WithIV(2)-WithoutIV(2),0,(WithIV_SD(2)+WithoutIV_SD(2))/2);
	[h(3) Ztest_pval(3,:)]=ztest(WithIV(3)-WithoutIV(3),0,(WithIV_SD(3)+WithoutIV_SD(3))/2);
	
	T=table(Variables,WithIV,WithIV_SD,WithoutIV,WithoutIV_SD,Ztest_pval)
	
	TEs=60e-3+tauGESSE(2:2:end);
	for k=1:length(Rs)
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		for j=1:length(TEs)
			[sigSE(j,k) tauSE sigSEev(j,k) sigSEiv(j,k)]=generate_signal(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','includeIV',true,'T2EV',Inf,'T2b0',Inf,'TE',TEs(j),'tau',0);
		end
	end
	
	figure;
	hold on;
	for k=1:length(Rs)
		plot(TEs.*1000,sigSEev(:,k),'-');
	end
	box;
	xlim([0 120]);
	ylim([0.8 1]);
	set(gca,'xtick',[-0:30:120]);
	set(gca,'ytick',[0.8:0.05:1]);
	legend('S_{EV}(R_c=1000\mum)','S_{EV}(R_c=50\mum)','S_{EV}(R_c=10\mum)','S_{EV}(R_c=5\mum)','location','south')
	grid on;
	title('Extravascular SE signal decay');
	xlabel('Echo time, TE (ms)');
	ylabel('Signal fraction (arb.)');
	
	figure;
	hold on;
	for k=1:length(Rs)
		plot(TEs.*1000,sigGESSEev(2:2:end,k)./sigSEev(:,k),'-');
	end
	box;
	xlim([0 120]);
	ylim([0.8 1]);
	set(gca,'xtick',[-0:30:120]);
	set(gca,'ytick',[0.8:0.05:1]);
	legend('S_{EV}(R_c=1000\mum)','S_{EV}(R_c=50\mum)','S_{EV}(R_c=10\mum)','S_{EV}(R_c=5\mum)','location','south')
	grid on;
	title('Extravascular GESSE signal decay normalised by SE signal decay');
	xlabel('Echo time, TE (ms)');
	ylabel('Signal fraction (arb.)');
	
