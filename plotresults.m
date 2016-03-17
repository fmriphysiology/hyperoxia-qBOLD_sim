function sig_ase=plotresults(p,storedPhase,TE);

	%NEED SOME ERROR CHECKING FOR CHOICE OF TE

	t=(p.deltaTE:p.deltaTE:p.TE*2)';
	
	if nargin>2
		p.TE=TE;
	end
	
	%error checking
	if or(mod(round(p.TE*1000),(round(p.deltaTE*2000)))>0,p.TE<=0)
		disp(['TE must be in the range ' num2str(p.deltaTE*2000) 'ms to ' num2str(t(end)*1000) 'ms in steps of ' num2str(p.deltaTE*2000) 'ms!']);
		ase=[];
		return;
	end
	
	%generate an ASE weighted signal
	TEind=find(round(t.*1000)==round(p.TE*1000),1,'first');
	
	for k=1:TEind+1
		ASEPhase(k,:)=sum(storedPhase(1:k-1,:),1)-sum(storedPhase(k:TEind,:),1);
	end
	tau_ase=(p.TE:-p.deltaTE*2:-p.TE)';
	sig_ase=abs(sum(exp(-i.*ASEPhase),2)./p.N);
	
	%generate a GESSE weighted signal
	TE2ind=find(round(t.*1000)==round(p.TE*500),1,'first');
	mask=repmat([ones(TE2ind,1); -ones(size(storedPhase,1)-TE2ind,1)],1,p.N);
	GESSEPhase=cumsum(storedPhase.*mask,1);
	tau_gesse=t-p.TE;
	sig_gesse=abs(sum(exp(-i.*GESSEPhase),2)./p.N);
	
	%plot signal curves
	figure(100);
	hold on;
	plot(tau_gesse.*1000,sig_gesse,'o-');
	box on;
	
	figure(101);
	hold on;
	plot(tau_ase.*1000,sig_ase,'o-');
	box on;

return;