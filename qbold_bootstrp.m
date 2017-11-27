function [sigEV tau params]=qbold_bootstrp(p,storedPhase,varargin)

	q=inputParser;
	addParameter(q,'TE',p.TE,@isnumeric); %target echo time
	addParameter(q,'Y',p.Y(end),@isnumeric); %target blood oxygen saturation
	addParameter(q,'Hct',p.Hct,@isnumeric); %target haematocrit
	addParameter(q,'permeable',false,@islogical); %are vessels permeable
	addParameter(q,'includeIV',false,@islogical); %include intravascular signal for ASE
	addParameter(q,'display',true,@islogical); %display plot of results
	addParameter(q,'T2EV',Inf,@isnumeric); %extravascular T2, defaults to infinite i.e. no effect 
	addParameter(q,'Vf',p.vesselFraction(1),@isnumeric);
	addParameter(q,'Yoff',0.95,@isnumeric); %blood oxygenation saturation saturation to match tissue susceptibility
	addParameter(q,'seq','ASE'); %blood oxygenation saturation saturation to match tissue susceptibility
	addParameter(q,'tau',[]);
	addParameter(q,'tau_cutoff',15e-3);
	parse(q,varargin{:});
	r=q.Results;

	t=(p.deltaTE:p.deltaTE:p.TE*2)';

	%error checking
	if or(mod(round(r.TE*1000),(round(p.deltaTE*2000)))>0,p.TE<=0)
		disp(['TE must be in the range ' num2str(p.deltaTE*2000) 'ms to ' num2str(t(end)*1000) 'ms in steps of ' num2str(p.deltaTE*2000) 'ms!']);
		sig_ase=[];
		return;
	end
	
	%scale for different Y values
	if length(p.R)>1
		Yscale=1;
	else
		Yscale=(1-r.Y)./(1-p.Y(end));
	end

	if strcmp(r.seq,'ASE')
		%generate an ASE weighted signal
		if isempty(r.tau) %error checking required here
			tau=(r.TE:-p.deltaTE*2:-r.TE)';
		else
			tau=r.tau;
			tau=tau(:);
		end
		
		TEind=r.TE/p.deltaTE;
		SEind=round((r.TE-tau)./(2.*p.deltaTE),0);
	
		for k=1:length(tau)
			Phase(k,:)=sum(storedPhase(1:SEind(k),:),1)-sum(storedPhase(SEind(k)+1:TEind,:),1);
		end	

	elseif strcmp(r.seq,'GESSE')
		%generate a GESSE weighted signal
		if isempty(r.tau)
			tau=(-r.TE./2:p.deltaTE:(p.TE.*2-r.TE))';
		else
			tau=r.tau;
		end
		TEind=round((r.TE+tau)./p.deltaTE,0);
		SEind=r.TE./(2*p.deltaTE);
	
		for k=1:length(tau)
			Phase(k,:)=sum(storedPhase(1:SEind,:),1)-sum(storedPhase(SEind+1:TEind(k),:),1);
		end		
	else
		disp('Sequence unknown');
		return;
	end

	numStepsInVessel=p.numStepsInVessel;
	
	switch(r.permeable)
		case true
			protonIndex=(1:5000);
		case false
			protonIndex=find(numStepsInVessel==0,5000,'first');
		otherwise
			protonIndex=find(numStepsInVessel==0,5000,'first');
	end
	
	if length(protonIndex)<5000
		sigEV=repmat(NaN,length(tau),1000);
		params=repmat(NaN,3,1000);
	else
		reverseStr = '';
		for k=1:1000
			sample_ind=round(mod(rand(5000,1).*1e6,4999),0)+1;
			sigEV(:,k)=calc_sig(Phase(:,protonIndex(sample_ind)),p,r);
			params(:,k)=calc_qbold_params(Phase(:,protonIndex(sample_ind)),p,r,tau)';
			%fprintf('.');
			
			percentDone = 100*k/1000;
			msg = sprintf('Percent done: %3.0f', percentDone); 
			fprintf([reverseStr, msg]);
			reverseStr = repmat(sprintf('\b'), 1, length(msg));
		end
		fprintf('\n');
	end
		
	function sigEV=calc_sig(Phase,p,r)

		Suscscale=(p.deltaChi0.*r.Hct.*(r.Yoff-r.Y))./(p.deltaChi0.*p.Hct.*(1-p.Y));
		sigEV=abs(sum(exp(-i.*Phase.*Suscscale),2)./length(Phase));
		
		if length(p.R)==1
			shapeEV=-log(sigEV)./p.vesselFraction;
			sigEV=exp(-r.Vf.*shapeEV);
		end
		
	return;

	function params=calc_qbold_params(Phase,p,r,tau)
	
		sigEV=calc_sig(Phase,p,r);
		
		tau_ind=find(tau>r.tau_cutoff);
		
		if length(tau_ind)<2
		
			r2pmeas=NaN;
			vmeas=NaN;
			e0meas=NaN;
		
		else
		
			X=[ones(size(tau_ind)) -tau(tau_ind)];
			a=X\log(sigEV(tau_ind));
		
			r2pmeas=a(2);
			vmeas=a(1)-log(sigEV(find(tau==0)));
			e0meas=(r2pmeas./vmeas)./(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*p.Hct);
		
		end
	
		params=[vmeas r2pmeas e0meas];
		
	return;
		
		