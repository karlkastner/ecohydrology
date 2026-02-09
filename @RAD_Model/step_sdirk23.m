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
function [z,stat] = step_sdirk23(obj,t,z0,dt)
	% solve nonlinear optimization problem
	g = (3+sqrt(3))/6;

	% first stage, intermediate value at t = g*dt
	dt_ = g*dt;
	obj.aux.rfun = @(z) z - z0 - dt_*obj.dz_dt(t+dt_,z);
	obj.aux.jfun = @jfun;
	[z1,stat] = obj.gauss_newton(z0);
	% check convergence
	if (stat.flag(end))
		z = z1;	
	else	
	
	% second stage, intermediate value at t = (1-g)*dt
	dz1dt = obj.dz_dt(t+dt_,z1);
	res0 = z0 + (1-2*g)*dt*dz1dt;
	dt_ = (1-g)*dt;
	obj.aux.rfun = @(z) z - g*dt*obj.dz_dt(t+dt_,z) - res0;

	% second stage, intermediate value at (1-g)*dt
	[z2,stat2] = obj.gauss_newton(z0);
	stat.flag  = [stat.flag, stat2.flag];
	stat.iter  = [stat.iter; stat2.iter];
	stat.linear_solver    = [stat.linear_solver, stat2.linear_solver];

	% step
	z   = z0 + 0.5*dt*(dz1dt + obj.dz_dt(t+(1-g)*dt,z2));
	end

	if (nargout()>1)
	% estimate the error
	tol = obj.time_integration_tolerance(z);
	
	% difference to midpoint scheme
	zmid = 0.5*(z + z0);
	dz   = (z - z0 - dt*obj.dz_dt(t+0.5*dt,zmid));	
	[dzmax,idmax] = max(abs(dz));
	% constant choosen so that for the linear ode dy/dt = -a*y the step
	% size (without safety margin) is chosen so that the error is equal to
	% the tolerance
	emax   = 0.59*dzmax;
	% safety margin is factored in later
	dt_opt = dt*cbrt(tol/emax);

	stat.rmse   = NaN;
	stat.maxe   = emax;
	stat.idmax  = idmax;
	stat.dt_opt = dt_opt;
	end

	function A = jfun(z)
		A = (obj.aux.I - g*dt*obj.jacobian(t+dt_,z,true));
	end % jfun
end % step_sdirk23

