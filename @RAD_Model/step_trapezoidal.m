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
function [z1,stat] = step_trapezoidal(obj,t,z0,dt)
	obj.aux.jfun = @jfun;
	obj.aux.rfun = @rfun;

	dzdt_0 = obj.dz_dt(t,z0);
	res0   = z0 + (0.5*dt)*dzdt_0;

	% solve nonlinear optimization problem
	[z1,stat] = obj.gauss_newton(z0);

	if (nargout()>1)
	% estimate the error
	tol = obj.time_integration_tolerance(z1);
		
	% as constant rates are linear, they do not contribute to the error
	% this avoids computing flux dependent constant rates
	dzdt_1 = obj.dz_dt(t+dt,z1);
	delta_dzdt     = (dzdt_1 - dzdt_0);
	[delta_dzdtmax,idmax] = max(abs(delta_dzdt));
	dzmax = dt*delta_dzdtmax;
	% this constant ensures that the estimated error is exact for the
	% linear ode dy/dt = a
	C      = 0.12;
	emax   = C*dzmax;
	% safety margin is factored in later
	dt_opt = dt*sqrt(tol/emax);

	stat.rmse   = NaN;
	stat.maxe  = emax;
	stat.idmax = idmax;
	stat.dt_opt = dt_opt;
	end

	function r = rfun(z)
		dz_dt = obj.dz_dt(t+dt,z);
		% (z1 - z0)/dt = q*A(z1) + (1-q)*A(z0)
		%r = z1 - z0 - ((1-q)*dt)*dz0_dt - dt*(q*dz1_dt + rhs);
		r = z - res0 - (0.5*dt)*dz_dt;
	end % rfun

	function A = jfun(z)
		Jt = obj.jacobian(t+dt,z,true);
		%if (obj.aux.surface_flow)
			A = (obj.aux.I - (0.5*dt)*Jt);
		%else
		%	A = (obj.aux.I_wo_surface_flow - (0.5*dt)*Jt);
		%end
	end % jfun
end % step_trapezoidal

