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
function [t, zz, out] = solve(obj)

	if (isa(obj.opt.solver,'function_handle'))
		solver_str = func2str(obj.opt.solver);
	else
		solver_str = obj.opt.solver;
	end

	switch (solver_str)
	case {'solve_split'}
		[t,zz,out] = obj.solve_split();
	case {'euler-forward'}
		[t,y,out]=obj.solve_euler_forward();
	otherwise % use matlab build-in solver
		out = struct('runtime',[]);
		dto = min(obj.opt.dto,obj.T);
		T = [0:dto:obj.T];
		tic();
		obj.init_advection_diffusion_matrix();
		out.runtime(1) = toc();
		tic();
		z0 = obj.opt.compute_class(obj.z0);
		[t,zz] = obj.opt.solver(@obj.dz_dt, T, z0, obj.odeopt);
		out.runtime(2) = toc();
		out.y_final = zz(end,:);
	end
end

