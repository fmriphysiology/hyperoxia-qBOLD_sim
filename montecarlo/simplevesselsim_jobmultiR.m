function simplevesselsim_jobmultiR(R1,R2)

	t=cputime;

	p=gentemplate;

	p.R=[R1 R2].*1e-6;
	p.Y=[0.6 0.6];
	p.Hct=[0.4 0.4];
	p.N=10000;
	p.universeScale=sqrt(25000);
	p.D=1e-9;
	p.vesselFraction=[0.5 0.5].*0.05;

	[spp p]=simplevesselsim(p);
		
	save(['../multi_vessel_radius_speedup/simvessim_res' num2str(R1) '-' num2str(R2) '.mat']);
	
	e=cputime-t;
	
	disp(['CPUtime (mins): ' num2str(e/60)]);

return;
