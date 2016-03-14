function plotresults(p,storedPhase);

	TE=p.TE.*1e3;
	nt=size(storedPhase,1);
	for k=1:nt/2
		ASEPhase(k,:)=sum(storedPhase(1:k,:),1)-sum(storedPhase(k+1:nt/2,:),1);
	end	
	tau_ase=(TE-4:-4:-TE)';
	ind=find((tau_ase>-34).*(tau_ase<61));
	
	storedPhase2=storedPhase.*repmat([ones(nt/4,1); -ones(3*nt/4,1)],1,p.N);
	GESSEPhase=cumsum(storedPhase2,1);
	tau_gesse=(2:2:TE*2)-TE;

	figure(100);
	hold on;
	plot(tau_gesse,(abs(sum(exp(-i.*GESSEPhase),2)./p.N)),'o-');
	grid on;
	box on;
	
	ind0=find(p.numStepsInVessel==0);
	figure(101);
	hold on;
	plot(tau_ase(ind),(abs(sum(exp(-i.*ASEPhase(ind,:)),2)./p.N)),'o-');
	%plot(tau_ase(ind),(abs(sum(exp(-i.*ASEPhase(ind,ind0)),2)./length(ind0))),'o-');
	grid on;
	box on;
	
	for k=1:100
		sig(:,k)=(abs(sum(exp(-i.*ASEPhase(ind,1:k*100)),2)./k*100));
	end
	sig_ss=sum((sig-repmat(sig(:,end),1,100)).^2);
	figure(102);
	hold on;
	plot((1:100).*100,sig_ss);
	grid on;
	box on;

%keyboard;

return;