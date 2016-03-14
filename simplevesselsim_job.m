function simplevesselsim_job(R)

	t=cputime;

	p=gentemplate;

	p.R=R*1e-6;
	p.N=10;

	[spp p]=simplevesselsim(p);
		
	save(['./svs_results' num2str(R) '.mat']);
	
	e=cputime-t;
	
	disp(['CPUtime (mins): ' num2str(e/60)]);

return;