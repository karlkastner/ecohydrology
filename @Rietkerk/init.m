% Mon 31 May 20:20:46 CEST 2021
function y0 = init(obj)
	dx      = obj.dx;
	%obj.dt  = 0.5*dx/rk.pmu.Dh;
	% initalize random number generator
	rng(obj.rng);

	% initial condition
	switch (obj.initial_condition)
	case {'random'}
			y0 = obj.random_state();
			% [b0,w0,h0] = obj.homogeneous_state(obj.p,state);
			% z0 = flat(repmat([b0,w0,h0],prod(obj.n),1));
	case {'colonized'}
			[b0,h0,w0] = obj.homogeneous_state([],0);
			o = ones(obj.n,1);
			y0 = [b0.*o; h0.*o; w0.*o];
			%fdx = rk.x <= meta.cL;
			% TODO no magic numbers
			y0(2) = 10;
	otherwise
		error('Rietkerk:init');
	end % swtich initial_condition

	% prepare spatially varying parameters
	f_C = fieldnames(obj.pmu);
	for jdx=1:length(f_C)
		field = f_C{jdx};
		sd    = obj.pss.(field);
		% scale
		% scaling standard deviations:
		% dx : devide   as perturbations are average in space when cells are larger
		% dt : multiply as perturbations should stay same over same time span,
		%      irrespective of number of time steps
		sd    = sd*sqrt(1/dx);
		%oname = [oname,'-sd_',field,'-',num2str(sd)];
		if (sd>0)
			obj.p.(field) = obj.pmu.(field).*gamrnd(1/sd^2,sd^2,obj.n,1);
		else
			obj.p.(field) = obj.pmu.(field);
		end
	%	obj.p.(field) = p.(field).mu*gr;
		%rw_gr = gamrnd(1/sda^2,sda^2,length(x),1);
		%Dh_gr = gamrnd(1/sda^2,sda^2,length(x),1);
	%	gr = meanfilt1(gr,nfa);
	%	for jdx=1:nfa
	%		gr = lowpass1d_implicit(gr,rhoa);
	%	end
	%	for jdx=1:nfh
	%		gr = highpass1d_implicit(gr,rhoa);
	%	end
	%	if (sda>0)
	%	gr=(gr-mean(gr));
	%	gr = 1+sda*gr/std(gr);
	%	end
	%	gr = gr(nfa/2+1:end);
	%	%gr = gr/mean(gr);
	%	p.(field)=p.(field)*gr;
	%	p.rw = p.rw*rw_gr;
	%	p.Dh = p.Dh*Dh_gr;
	end
	%rk.p = p;

%	obj.x = linspace(0,obj.n(1),obj.L(1));
	switch (length(obj.n))
	case {1}
		if (obj.pmu.vh ~= 0)
			obj.D1  = derivative_matrix_1_1d(obj.n,obj.L,-sign(obj.p.vh(1)),'c');
		else
			obj.D1 = 0;
		end
		obj.D1c = derivative_matrix_1_1d(obj.n,obj.L,2,'c');
		obj.D2  = derivative_matrix_2_1d(obj.n,obj.L,obj.order,'circular');
	case {2}
%		obj.y   = linspace(0,obj.n(2),obj.L(2));
		obj.D1  = [];
		obj.D1c = [];
		[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(obj.n,obj.L,obj.order,'circular');
		obj.D2 = D2x+D2y;
	otherwise
		error('n and L must have both 1 or 2 elements')
	end
	n = prod(obj.n);
	obj.Z = sparse(n,n);
	obj.I = speye(n);

end

