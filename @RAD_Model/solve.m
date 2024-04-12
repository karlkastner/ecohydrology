% 2021-07-06 10:56:59.150856622 +0200
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
function [to, zo, out] = solve(obj)
	% convert solver to string for switch-case
	if (isa(obj.opt.solver,'function_handle'))
		solver_str = func2str(obj.opt.solver);
	else
		solver_str = obj.opt.solver;
	end

	%
	% call solver
	%
	switch (solver_str)
	case {'euler_forward','solve_implicit','solve_split', ...
		'step_react_advect_diffuse_erk', ...
		'step_integrating_factor'}
		[to,zo,out] = obj.solve_step();
	otherwise % use matlab build-in solver
		dto = min(obj.opt.dto,obj.T);
		T = [0:dto:obj.T];
		timer=tic();
		z0 = obj.opt.compute_class(obj.z0);
		[to,zo] = obj.opt.solver(@obj.dz_dt, T, z0, obj.odeopt);
		obj.out.runtime(2) = toc(timer);
		obj.out.y_final = zz(end,:);
	end % switch solver_str
	out = obj.out;
end % solve

