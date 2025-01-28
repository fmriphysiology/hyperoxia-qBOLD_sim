function [params paramssd]=calc_qbold_params(p,sig,tau,tau_cutoff,Hct)

	if nargin<4
		tau_cutoff=15e-3;
	end
	
	if nargin<5
		Hct=p.Hct;
	end
	
	tau=tau(:);
	tau_ind=find(tau>=tau_cutoff);
	
	if length(tau_ind)<2
	
		r2pmeas=NaN;
		vmeas=NaN;
		e0meas=NaN;
		s0meas=NaN;
		
		r2pmeassd=NaN;
		vmeassd=NaN;
		e0meassd=NaN;
		s0meassd=NaN;
		
	else
	
		X=[ones(size(tau_ind)) -tau(tau_ind) ones(size(tau_ind))];
		X=[1 0 0; X];
		Y=sig(tau_ind);
		Y=[sig(tau==0); Y];
		[a asd mse s]=lscov(X,log(Y));
	
		r2pmeas=a(2);
		vmeas=a(3);
		e0meas=(r2pmeas./vmeas)./(4/3*pi*p.gamma.*p.B0.*p.deltaChi0.*Hct);
		s0meas=a(1);
		
		r2pmeassd=asd(2);
		vmeassd=asd(3);
		e0meassd=sqrt((r2pmeassd./r2pmeas).^2+(vmeassd./vmeas).^2+2.*s(1,2)./(r2pmeas.*vmeas))*e0meas;
		s0meassd=asd(1);
	
	end

	params=[r2pmeas vmeas e0meas s0meas]';
	paramssd=[r2pmeassd vmeassd e0meassd s0meassd]';
	