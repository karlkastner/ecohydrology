% Mon 31 May 20:20:46 CEST 2021
% Karl Kästner, Berlin
%
%% initialize all variables
%
function y0 = init(obj)
	dx      = obj.dx;
	obj.fx  = fourier_axis(obj.x);
	if (obj.ndim>1)
	obj.fy  = fourier_axis(obj.x(2));
	end
	% initalize random number generator
	rng(obj.rng);

	% initial condition
	if (isnumeric(obj.initial_condition))
		y0 = obj.initial_condition;
	else
	switch (obj.initial_condition)
	case {'random'}
			y0 = obj.random_state();
			% [b0,w0,h0] = obj.homogeneous_state(obj.p,state);
			% z0 = flat(repmat([b0,w0,h0],prod(obj.n),1));
	case {'colonize'}
			[b0,h0,w0] = obj.homogeneous_state([],0);
			o = ones(obj.n,1);
			y0 = [b0.*o; h0.*o; w0.*o];
			% TODO no magic numbers
			y0(2) = 10;
	otherwise
		error('Rietkerk:init');
	end % swtich initial_condition
	end

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
		n = obj.n;
		if (length(n)<2)
			n(2) = 1;
		end

		if (sd>0)
			if (~obj.opt.gbm)
				sd_averaged   = sd*sqrt(1/prod(dx)); %dx.^obj.ndim);
				obj.p.(field) = flat(obj.pmu.(field).*gamrnd(1/sd_averaged^2,sd_averaged^2,n));
			else
				% by definition, the moments of the (g)-bm-(bridge)
				% do not depend on dx but only on L, so dx does
				% not have to be rescaled
				%
				% TODO this is for 1d only
				% TODO simply make log(mu_a) = bar log(x)
				%		   
				% determine parameter
				[gmu,gsd] = gbm_moment2par(1,sd,obj.x(end)-obj.x(1));
				% simulate geometric brownian bridge
				z = gbm_bridge(obj.x,gmu,gsd,1,1);
				obj.p.(field) = flat(obj.pmu.(field))*z;
			end
		else
			obj.p.(field) = obj.pmu.(field);
		end
	end

	% expand scalar advection and diffusion coefficients (isotropic)
	if (obj.ndim > 1)
		f_C = {'eh'}; %,'eb','ew'};
		for idx=1:length(f_C)
			field = f_C{idx};
			if (1 == length(obj.p.(field)))
				obj.p.(field)(2) = obj.p.(field)(1);
			end
		end
		% the velocity of the second axis is set to zero here
		if (1 == length(obj.p.vh))
			obj.p.vh(2) = 0;
		end
	end

	% prepare discretization matrices
	switch (length(obj.n))
	case {1}
		if (obj.pmu.vh ~= 0)
			obj.D1x  = derivative_matrix_1_1d(obj.n,obj.L,-sign(obj.p.vh(1)),obj.bc{1});
		else
			obj.D1x = spzeros(obj.n(1));
		end
		obj.D1xc = derivative_matrix_1_1d(obj.n,obj.L,2,obj.bc{1});
		obj.D2x  = derivative_matrix_2_1d(obj.n,obj.L,obj.order,obj.bc{1});
		obj.D1   = obj.D1x;
		obj.D1c  = obj.D1xc;
		obj.D2   = obj.D2x;
	case {2}
		if (obj.p.vh(1) ~= 0)
			D1x  = derivative_matrix_1_1d(obj.n(1),obj.L(1),-sign(obj.p.vh(1)),obj.bc{2});
			%obj.D1x = kron(D1x,speye(obj.n(2)));
			obj.D1x = kron(speye(obj.n(2)),D1x);
		else
			obj.D1x = spzeros(prod(obj.n));
			%parse(prod(obj.n));
		end
		if (obj.p.vh(2) ~= 0)
			D1y  = derivative_matrix_1_1d(obj.n(2),obj.L(2),-sign(obj.p.vh(2)),obj.bc{2});
			%obj.D1y = kron(speye(obj.n(1)),D1y);
			obj.D1y = kron(D1y,speye(obj.n(1)));
		else
			%obj.D1y = sparse(prod(obj.n));
			obj.D1y = spzeros(prod(obj.n));
		end
		obj.D1xc = 0;
		obj.D1yc = 0;
		[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(obj.n,obj.L,obj.order,obj.bc);
		obj.D2  = D2x+D2y;
		obj.D2x = D2x;
		obj.D2y = D2y;
		%obj.D1  = obj.D1x+obj.D1y;
	otherwise
		error('n and L must have both 1 or 2 elements')
	end
	n     = prod(obj.n);
	obj.Z = sparse(n,n);
	obj.I = speye(n);

end

