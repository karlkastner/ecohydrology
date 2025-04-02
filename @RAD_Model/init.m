% Mon 31 May 20:20:46 CEST 2021
% Karl Kästner, Berlin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% function obj = init(obj)
function obj = init(obj)
	n = obj.nx;
	L = obj.L;

	% initalize random number generator
	if (~isempty(obj.opt.rng))
		rng(obj.opt.rng);
	end

	% spatially varying parameters
	field_C = fieldnames(obj.pmu);
	for idx=1:length(field_C)
		% check if the parameters has not yet been defined
		% as constructor argument
		if (~isfield(obj.p,field_C{idx}) || isempty(obj.p.(field_C{idx})))
			mu = obj.pmu.(field_C{idx});
			sd = obj.pss.(field_C{idx});
			sl = obj.psl.(field_C{idx});
			dist = obj.psdist.(field_C{idx});
			%n = prod(obj.nx);
			[obj.p.(field_C{idx}),C,S] = generate(dist,mu,sd,sl,obj.nx); %,field_C{idx});
			if (~isempty(C))
				obj.psC.(field_C{idx}) = C;
			end
			if (~isempty(S))
				obj.psS.(field_C{idx}) = S;
			end
		end
	end

	% initial condition
	if (isstruct(obj.initial_condition))
		z0 = [];
		for idx=1:obj.nvar
			dist = obj.initial_condition.dist{idx};
			mu = obj.initial_condition.mu(idx);
			ic = obj.initial_condition
			% random
			if (~isempty(obj.initial_condition.sd))
				sd = obj.initial_condition.sd(idx);
				sl = obj.initial_condition.sl(idx);
				z0i = flat(generate(dist,mu,sd,sl,obj.nx));
				if (isscalar(z0i))
					z0i = z0i*ones(prod(obj.nx),1);
				end
			elseif (isfield(obj.initial_condition,'fc') && ~isempty(obj.initial_condition.fc))
				fc = obj.initial_condition.fc
				c = obj.initial_condition.c(idx)
				s = obj.initial_condition.s(idx)
				% TODO 2D
				z0i = mu + s*sin(2*pi*fc*cvec(obj.x)) + c*cos(2*pi*fc*cvec(obj.x));
			end
			z0 = [z0;z0i];
		end
		obj.z0 = z0;
	elseif (isnumeric(obj.initial_condition))
		obj.z0 = obj.initial_condition;
	else
		obj.z0 = eval(obj.initial_condition);
	end

	% initialize_solver
	obj.init_solve();

function [x,C,S] = generate(dist,mu,sd,sl,n)
		C = [];
		S = [];
		dx = obj.dx;
		% scaling standard deviations:
		% dx : devide   as perturbations are average in space when cells are larger
		% dt : multiply as perturbations should stay same over same time span,
		%      irrespective of number of time steps
		sd = sd./sqrt(prod(dx));
		if (isempty(dist)||0==sd)
		if (0==sd)
			x = mu;
		else
			error('when variance is not zero, a distribution has to be specified');
		end
		else
		switch (dist)
		case {'uniform'}
			[a,b] = uniform_moment2par(mu,sd);
			x = uniformrnd(a,b,prod(n),1);
		case {'normal'}
			x = normrnd(mu,sd,prod(n),1);
		case {'gamma'}
			[a,b] = gampdf_moment2par(mu,sd);
			x = gamrnd(a,b,prod(n),1);
			% x = flat(mu.*gamrnd(1/sd^2,sd^2,n);
		case {'lognormal'}
			[a,b] = lognpdf_moment2par(mu,sd);
			x = lognrnd(a,b,prod(n),1);
		case {'exp'}
			x = exprnd(mu,prod(n),1);
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
					error('gbb not yet implemented in 2d');
				end
		case {'geometric-ornstein-uhlenbeck'}
			if (1 == obj.ndim)
				error('1d not yet implemented');
			else
				% determine parameter
				%[lmu,lmu_,lsd,lsd_,ltheta] = lognpdf_moment2par_correlated(mu,mu,sd,sd,sl);
				%ltheta = geometric_ou_correlation_length_of_z(sl,ltheta);
				oi = obj.opt.heterogeneity_integration_order;
				oL = obj.opt.heterogeneity_oversampling_factor_spectral;
				pw = obj.opt.heterogeneity_p_window;
				% simulate an instantiation of the process
				theta = sl;
				[x,C,S] = geometric_ou_2d_grid_cell_averaged_generate(mu,sd,theta,obj.L,obj.nx,oi,oL,pw);
				% matrix to vector
				x  = flat(x);
			end
		case {'geometric-pink'}
			if (1 == obj.ndim)
				error('1d not yet implemented');
			else
				% pink noise with unit variance
				e = pink_noise_2d(n,L);
				% rescale to desired mean and variance
				[lmu,lsd] = lognpdf_moment2par(mu,sd);
				% transform to log-normal
				e = exp(lmu+lsd*e);
				x = flat(e);
			end
		otherwise
			error('unimplemented distribution');
		end
		end
end % generate

end %

