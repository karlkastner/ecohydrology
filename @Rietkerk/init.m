% Mon 31 May 20:20:46 CEST 2021
% Karl Kästner, Berlin
%
%% initialize all variables
%
function init(obj)

	% initalize random number generator
	if (~isempty(obj.opt.rng))
		rng(obj.opt.rng);
	end

	% initial condition
	if (isnumeric(obj.opt.initial_condition))
		z0 = obj.opt.initial_condition;
	else
	switch (obj.opt.initial_condition)
	case {'random'}
			z0 = obj.random_state();
			% [b0,w0,h0] = obj.homogeneous_state(obj.p,state);
			% z0 = flat(repmat([b0,w0,h0],prod(obj.n),1));
	case {'colonize'}
			[b0,h0,w0] = obj.homogeneous_state([],0);
			o = ones(obj.n,1);
			z0 = [b0.*o; h0.*o; w0.*o];
			% TODO no magic numbers
			z0(2) = 10;
	otherwise
		error('Rietkerk:init');
	end % swtich initial_condition
	end

	% prepare spatially varying parameters
	f_C = fieldnames(obj.pmu);
	for jdx=1:length(f_C)
		field = f_C{jdx};
		% check if it has been externally defined already
		if (~isfield(obj.p,field) || isempty(obj.p.(field)))
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
			switch (obj.opt.heterogeneity_model)
			case {'gamma'}
				sd_averaged   = sd*sqrt(1./prod(obj.dx)); %dx.^obj.ndim);
				obj.p.(field) = flat(obj.pmu.(field).*gamrnd(1/sd_averaged^2,sd_averaged^2,n));
			case {'geometric-brownian-bridge'}
				if (1 == obj.ndim)
				% by definition, the moments of the (g)-bm-(bridge)
				% do not depend on dx but only on L, so dx does
				% not have to be rescaled
				%		   
				% determine parameter
				[gmu,gsd] = gbm_moment2par(1,sd,obj.x(end)-obj.x(1));
				% simulate geometric brownian bridge
				z = gbm_bridge(obj.x,gmu,gsd,1,1);
				obj.p.(field) = obj.pmu.(field)*z;
				else
					error('not yet implemented');
				end
			case {'geometric-ornstein-uhlenbeck'}
				if (1 == obj.ndim)
					error('not yet implemented');
				else
					% determine parameter
					[lmu,lsd] = logn_moment2par(obj.pmu.(field),obj.pss.(field));
					theta = obj.pssl.(field);
					% TODO no magic numbers
					m = obj.opt.heterogeneity_integration_order;
					z = geometric_ar1_2d_grid_cell_averaged_generate(lmu,lsd,theta,obj.L,obj.n,m);
					obj.p.(field) = z;
				end
			otherwise
				error('not yet implemented');
			end
		else
			obj.p.(field) = obj.pmu.(field);
		end
		else
			% use prespecified value and do nothing
			% TODO check size
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
	% TODO only when required by solver
	switch (length(obj.n))
	case {1}
		if (obj.pmu.vh ~= 0)
			obj.aux.D1x  = derivative_matrix_1_1d(obj.n,obj.L,-sign(obj.p.vh(1)),obj.bc{1});
		else
			obj.aux.D1x = spzeros(obj.n(1));
		end
		obj.aux.D1xc = derivative_matrix_1_1d(obj.n,obj.L,2,obj.bc{1});
		obj.aux.D2x  = derivative_matrix_2_1d(obj.n,obj.L,obj.opt.order,obj.bc{1});
		obj.aux.D1   = obj.aux.D1x;
		obj.aux.D1c  = obj.aux.D1xc;
		obj.aux.D2   = obj.aux.D2x;
	case {2}
		if (obj.p.vh(1) ~= 0)
			D1x  = derivative_matrix_1_1d(obj.n(1),obj.L(1),-sign(obj.p.vh(1)),obj.bc{2});
			obj.aux.D1x = kron(speye(obj.n(2)),D1x);
		else
			obj.aux.D1x = spzeros(prod(obj.n));
		end
		if (obj.p.vh(2) ~= 0)
			D1y  = derivative_matrix_1_1d(obj.n(2),obj.L(2),-sign(obj.p.vh(2)),obj.bc{2});
			obj.aux.D1y = kron(D1y,speye(obj.n(1)));
		else
			obj.aux.D1y = spzeros(prod(obj.n));
		end
		obj.aux.D1xc = 0;
		obj.aux.D1yc = 0;
		[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(obj.n,obj.L,obj.opt.order,obj.bc);
		obj.aux.D2  = D2x+D2y;
		obj.aux.D2x = D2x;
		obj.aux.D2y = D2y;
		%obj.D1  = obj.D1x+obj.D1y;
	otherwise
		error('n and L must have both 1 or 2 elements')
	end
	n     = prod(obj.n);
	obj.aux.Z = sparse(n,n);
	obj.aux.I = speye(n);

	obj.z0 = z0;

	% sanity check
	field_C = {'vh','eh','eb','ew'};
	for idx=1:length(field_C)
		if (length(obj.pmu.(field_C{idx})) ~= obj.ndim)
			error(['Dimension of parameter ',field_C{idx},' does not match']);
		end
	end
end

