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
%% two stage euler, with richardson extrapolation,
%% second order accurate
%%
function [z,stat] = step_euler2x(obj,t,z0,dt)
	% solve nonlinear optimization problem
	% first half step
	obj.aux.rfun = @(z1) rfun(t,0.5*dt,z0,z1);
	obj.aux.jfun = @(z1) jfun(t,0.5*dt,z1);
	[z1,stat] = obj.gauss_newton(z0);
	% second hal step
	obj.aux.rfun = @(z2) rfun(t+0.5*dt,0.5*dt,z1,z2);
	obj.aux.jfun = @(z2) jfun(t+0.5*dt,0.5*dt,z2);
	[z2,stat2] = obj.gauss_newton(z1);
	stat.flag  = [stat.flag, stat2.flag];
	stat.iter  = [stat.iter; stat2.iter];
	stat.linear_solver    = [stat.linear_solver, stat2.linear_solver];
	% coarse step, use z2 as initial condition
	obj.aux.rfun = @(z2_) rfun(t,dt,z0,z2_);
	obj.aux.jfun = @(z2_)  jfun(t,dt,z2_);
	[z2_,stat3] = obj.gauss_newton(z0);
	stat.flag  = [stat.flag, stat3.flag];
	stat.iter  = [stat.iter; stat3.iter];
	stat.linear_solver    = [stat.linear_solver, stat3.linear_solver];
	% error estimate
	p = 1;
	e = (z2_ - z2)./(2.^p-1);
	% richardson extrapolation
	z = z2 - e;

	if (nargout()>1)
	% estimate the error
	tol = obj.time_integration_tolerance(z);

	% c.f. johnson 1988
	% TODO allow for different error norms
	[emax,idmax] = max(abs(e));
	% the constant is that of the error when solving the
	% linear ode dz/dt = -a z
	% this is about twice as large as trapezoidal or midpoint
	C = 0.24;
	dt_opt = C*dt*sqrt(tol/emax);
	
	stat.rmse   = NaN;
	stat.maxe   = emax;
	stat.idmax  = idmax;
	stat.dt_opt = dt_opt;
	end

	function r = rfun(t,dt,z0,z1)
		r = z1 - z0 - dt*obj.dz_dt(t+dt,z1);
	end % rfun

	function A = jfun(t,dt,z)
		Jt = obj.jacobian(t+dt,z,true);
		A = (obj.aux.I - dt*Jt);
	end % jfun
end % step_euler_double

