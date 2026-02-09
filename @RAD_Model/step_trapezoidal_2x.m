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
function [z,stat] = step_trapezoidal_2x(obj,t,z0,dt)
	dz_dt0 = obj.dz_dt(t,z0);
	res0   = z0 + (0.25*dt)*dz_dt0;
	obj.aux.jfun = @(z) (obj.aux.I - (0.25*dt)*obj.jacobian(t+dt/2,z,true));
	obj.aux.rfun = @(z) (z - res0 - (0.25*dt)*obj.dz_dt(t+dt/2,z));
	% solve nonlinear optimization problem
	[z1,stat] = obj.gauss_newton(z0);
	dz_dt1 = obj.dz_dt(t+dt/2,z1);
	res1   = z1 + (0.25*dt)*dz_dt1;
	obj.aux.jfun = @(z) (obj.aux.I - (0.25*dt)*obj.jacobian(t+dt,z,true));
	obj.aux.rfun = @(z) (z - res1 - (0.25*dt)*obj.dz_dt(t+dt,z));
	% solve nonlinear optimization problem
	[z2,stat2] = obj.gauss_newton(z1);
	stat.flag  = [stat.flag, stat2.flag];
	stat.iter  = [stat.iter; stat2.iter];
	stat.linear_solver    = [stat.linear_solver, stat2.linear_solver];

	% TODO pass z2 as initial guess
	res0   = z0 + (0.5*dt)*dz_dt0;
	obj.aux.jfun = @(z) (obj.aux.I - (0.5*dt)*obj.jacobian(t+dt,z,true));
	obj.aux.rfun = @(z) (z - res0 - (0.5*dt)*obj.dz_dt(t+dt,z));
	[z2_,stat3] = obj.gauss_newton(z0);
	stat.flag  = [stat.flag, stat3.flag];
	stat.iter  = [stat.iter; stat3.iter];
	stat.linear_solver    = [stat.linear_solver, stat3.linear_solver];
	
	% error estimate
	p = 2;
	e = (z2_ - z2)./(2.^p-1);
	% richardson extrapolation
	z = z2 - e;

	if (nargout()>1)
	% estimate the error
	tol = obj.time_integration_tolerance(z);
		
	% as constant rates are linear, they do not contribute to the error
	% this avoids computing flux dependent constant rates
	%obj.aux.ignore_constant_rates = true;
	%dzdt_0 = obj.dz_dt(t,z0);
	%dzdt_1 = obj.dz_dt(t+dt,z1);
	%obj.aux.ignore_constant_rates = false;
	delta_dzdt     = e; %(dzdt_1 - dzdt_0);
	[delta_dzdtmax,idmax] = max(abs(delta_dzdt));
	dzmax = dt*delta_dzdtmax;
	% this constant ensures that the estimated error is exact for the
	% linear ode dy/dt = a
	% TODO this constant is for trap, not trap2x
	C      = 0.12;
	emax   = C*dzmax;
	% safety margin is factored in later
	dt_opt = dt*cbrt(tol/emax);

	stat.rmse   = NaN;
	stat.maxe  = emax;
	stat.idmax = idmax;
	stat.dt_opt = dt_opt;
	end


end % step_trapezoidal

