% Thu  7 Dec 16:41:17 CET 2023
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
function dt = adapt_time_step(obj,t,dt_old)
	% max time step satisfying the error constraint
	dt = obj.aux.stat.dt_opt;

	dt = dt.*obj.opt.time_integration.dt_safety_factor;

 	% limit decrease, this is usually unlimetted
	dt = max(dt,dt_old*obj.opt.time_integration.dt_min_scale_down);

	% limit increase
	dt = min(dt,dt_old*obj.opt.time_integration.dt_max_scale_up);

	% limit at lower bound, this is usually only limited to 0
	dt = max(dt,obj.opt.time_integration.dt_min);

	% limit at upper bound
	dt = min(dt,obj.opt.time_integration.dt_max);

	% limit by condition number
	% linsolver_tol = 1/c = 1/(1 + q*dt*max(eig(J))
	tol = obj.opt.linear_solver.tol;
	obj.aux.dt_max_prec = obj.opt.time_integration.dt_max_prec_scale*(1/tol-1)/(obj.aux.q*obj.aux.max_eig_J);
	dt  = min(dt,obj.aux.dt_max_prec);

	% avoid too many inner steps
	if (~isempty(obj.opt.nlsolver.name))
	switch (obj.opt.nlsolver.name)
	case {'gauss-newton'}
	switch (obj.opt.linear_solver.name)
	% TODO mg-java was removed
	case {'multigrid-java'}
		iter_max = obj.opt.inner2_maxiter_;
		if (~isempty(iter_max))
			% the last one is the hardest always
			iter = obj.aux.stat.linear_solver.iter(end);
			dt = min(dt,dt_old*iter_max/iter);
		end
	end % switch i2
	end % switch innersolver
	end
end

