function simplevesselsim_jobR(R)

	t=cputime;

	p=gentemplate;

	p.R=R*1e-6;
	p.N=10000;
	p.universeScale=sqrt(25000);
	p.D=1e-9;

	[spp p]=simplevesselsim(p);
		
	save(['../simbloodsim_res' num2str(R) '.mat']);
	
	e=cputime-t;
	
	disp(['CPUtime (mins): ' num2str(e/60)]);

return;
