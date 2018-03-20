function figure_investigate_te(simdir)

	Rs=[1:10 20:10:100 200:100:1000];
	
	%change DBV, fix OEF
	
	%ASE signal with D=1, Vf=3, Y=60
	Vf=0.03;
	Y=0.6;

	TE=80e-3;
	for k=1:length(Rs)
		disp(['Step ' num2str(k) ' of ' num2str(length(Rs))]);
		load([simdir 'single_vessel_radius_D1-0Vf3pc/simvessim_res' num2str(Rs(k))]);
		[sigASED1V3Y60(:,:,k) tauASE paramsASED1V3Y60(:,:,k)]=qbold_bootstrp(p,spp,'display',false,'Vf',Vf,'Y',Y,'seq','ASE','TE',TE);
	end
