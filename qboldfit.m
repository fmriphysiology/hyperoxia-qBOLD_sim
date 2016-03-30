function [vmeas r2pmeas e0meas ssq]=qboldfit(p,storedPhase,TE,Y)

	t=(p.deltaTE:p.deltaTE:p.TE*2)';
	deltaChi0=0.264e-6;
	Hct=0.4;
	
	if nargin>2
		if ~isempty(TE)
			p.TE=TE;
		end
	end
	
	if nargin<4
		Y=p.Y(end);
	end
		
	%error checking
	if or(mod(round(p.TE*1000),(round(p.deltaTE*2000)))>0,p.TE<=0)
		disp(['TE must be in the range ' num2str(p.deltaTE*2000) 'ms to ' num2str(t(end)*1000) 'ms in steps of ' num2str(p.deltaTE*2000) 'ms!']);
		ase=[];
		return;
	end

	%scale for different Y values
	Yscale=(1-Y)./(1-p.Y(end));
	
	%generate an ASE weighted signal
	TEind=find(round(t.*1000)==round(p.TE*1000),1,'first');
	
	for k=1:TEind+1
		ASEPhase(k,:)=sum(storedPhase(1:k-1,:),1)-sum(storedPhase(k:TEind,:),1);
	end
	tau_ase=(p.TE:-p.deltaTE*2:-p.TE)';
	sig_ase=abs(sum(exp(-i.*ASEPhase(:,find(p.numStepsInVessel==0)).*Yscale),2)./sum(p.numStepsInVessel==0));
	%sig_ase=abs(sum(exp(-i.*ASEPhase),2)./p.N);
	sum(p.numStepsInVessel==0)
	
	tau_ind=find(tau_ase>15e-3);
	X=[ones(size(tau_ind)) -tau_ase(tau_ind)];
	a=X\log(sig_ase(tau_ind));
	
	ssq=sum((log(sig_ase(tau_ind))-X*a).^2)./length(tau_ind);
	
	r2pmeas=a(2);
	vmeas=a(1)-log(sig_ase(find(tau_ase==0)));
	e0meas=(r2pmeas./vmeas)./(4/3*pi*p.gamma.*p.B0.*deltaChi0.*Hct);
	
	%keyboard;