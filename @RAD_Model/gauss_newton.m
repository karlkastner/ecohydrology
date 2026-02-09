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
% TODO merge with common gn-function
function [z,stat] = gauss_newton(obj,z)
	stat = struct('rmse',NaN ...
		      ,'dt0', NaN ...
		      ,'flag', 0 ...
		      ,'iter', NaN ...
		      ,'ls', struct() ...
		      ,'f0', NaN ...
		      ,'f', NaN  ...
		      ,'linear_solver', struct('flag',[],'relres',[],'iter',[]) ...
		      ,'flag2', NaN ...
		      ,'fs', NaN ...
		      ,'runtime', NaN ...
		      ,'dt', NaN ...
		      ,'maxe', NaN ...
		      ,'idmax', NaN ...
		      ,'dt_opt', NaN ...
			);


	% TODO use generic gn-function
	f    = [];
	k = 0;
	while (1)
		k = k+1;
		stat.iter = k;
		timeri = tic();
		[stat.f0(k),g,A] = gnfun(z);
		if (~isfinite(stat.f0(k)) && obj.opt.exception_on_nan)
			error('Objective function is not finite');
		end
		
		% TODO distinguish abstol from reltol
		% TODO convergence criterion based on gradiend
		%if (stat.f0(k) <= obj.opt.nlsolver_tol*obj.opt.nlsolver_tol) % || g'*g < 1e-12)
		if (g'*g < obj.opt.nlsolver.abstol*obj.opt.nlsolver.abstol)
			break;
		end
		if (k > obj.opt.nlsolver.maxiter)
			warning('nlsolver reached maxiter')
			stat.flag(k) = -1;
			return;
		end
		[dz,stat.linear_solver(k)] = solve_linear_system(A,g,obj.opt.linear_solver);
		if (stat.linear_solver(k).flag)
			stat.flag(k) = -2;
			return;
		end
		dz = -dz;
		a_ = 1;
		%dz = max(dz,-z);
		[as,stat.ls.flag(k),stat.fs(k),stat.ls.iter(k)] = line_search_backtracking(z,stat.f0(k), g, @gnfun,dz,a_,obj.opt.nlsolver.line_search_maxiter);
		z = z + as*dz;
		%if (obj.opt.nlsolver.line_search_ispos && min(z)<0)
		% so interestingly, the value can be negative after the first GN step, but is positive after the second
		if (obj.opt.nlsolver.line_search_ispos && min(z)<-obj.opt.linear_solver.tol)
			stat.flag(k) = -4;
			warning('negative value');
			break;
		end
		stat.ls.as(k) = as;
		stat.runtime(k) = toc(timeri);
		if (stat.ls.flag(k))
			warning('Line search did not converge');
			stat.flag(k) = -3;
			break;
		end
	end % while no convergence
	%otherwise
	%	error('unimplemented solver')
	%	end % switch opt.nlsolver, case newton

function [f,g,A] = gnfun(z)
	% residual r is equal to gradient g
	g = obj.aux.rfun(z);

	% objective function value
	f = 0.5*(g'*g);

	if (nargout()>2)
		A = obj.aux.jfun(z);
	end
end % gnfun

end % gauss_newton

