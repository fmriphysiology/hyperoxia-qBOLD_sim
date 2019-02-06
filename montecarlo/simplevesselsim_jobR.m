function simplevesselsim_jobR(R)

	t=cputime;

	p=gentemplate;

	p.R=R*1e-6;
	p.N=10000;
	p.universeScale=sqrt(25000);
	p.D=1e-9;
	p.vesselFraction=0.03;
	p.Y=0.6;

	[spp p]=simplevesselsim(p);
		
	save(['../single_vessel_radius_D1-0Vf3pc_dist/simvessim_res' num2str(R) '.mat']);
	
	e=cputime-t;
	
	disp(['CPUtime (mins): ' num2str(e/60)]);

return;
