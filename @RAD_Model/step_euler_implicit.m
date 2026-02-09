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
%
function [z1,stat] = step_euler_implicit(obj,t,z0,dt)
	if (true)
	obj.aux.jfun = @jfun;
	obj.aux.rfun = @rfun;
	% solve nonlinear optimization problem
	[z1,stat] = obj.gauss_newton(z0);
	else
		% TODO inhomogeneous boundary
		obj.aux.rhs  = z0;
		obj.aux.afun = @(t,h) (obj.aux.I - dt*obj.fixed_point_matrix(t,h));
		[z1,stat]=obj.fixed_point_iteration(z0);
	end

	if (nargout()>1)

	% estimate the error
	tol = obj.time_integration_tolerance(z1);

	% c.f. johnson 1988
	dz   = (z1 - z0);
	% TODO allow for different error norms
	[emax,idmax] = max(abs(dz));
	% the constant is that of the truncation error when solving the
	% linear ode dz/dt = -a z
	C      = 0.072;
	emax   = C*emax; 
	dt_opt = dt*tol/emax;
	
	stat.rmse   = NaN;
	stat.maxe   = emax;
	stat.idmax  = idmax;
	stat.dt_opt = dt_opt;
	end

	function r = rfun(z1)
		r = z1 - z0 - dt*obj.dz_dt(t+dt,z1);
	end % rfun

	function A = jfun(z1)
		Jt = obj.jacobian(t+dt,z1,true);
		%if (obj.aux.surface_flow)
			A = (obj.aux.I - dt*Jt);
		%else
		%	A = (obj.aux.I_wo_surface_flow - dt*Jt);
		%end
	end % jfun
end % step_euler_implicit

