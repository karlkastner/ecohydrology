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
function [z,stat] = step_midpoint(obj,t,z0,dt)
	obj.aux.jfun = @jfun;
	obj.aux.rfun = @rfun;

	% solve nonlinear optimization problem
	[z,stat] = obj.gauss_newton(z0);

	if (nargout()>1)

	% estimate the error
	tol = obj.time_integration_tolerance(z);

	% as constant rates are linear, they do not contribute to the error
	% this avoids computing flux dependent constant rates
	obj.aux.ignore_constant_rates = true;
	dzdt_0 = obj.dz_dt(t,z0);
	dzdt_1 = obj.dz_dt(t,z);
	obj.aux.ignore_constant_rates = false;
	delta_dzdt  = (dzdt_1 - dzdt_0);
	[delta_dzdtmax,idmax] = max(abs(delta_dzdt));
	dzmax = dt*delta_dzdtmax;
	% this constant ensures that the estimated error is exact for the
	% linear ode dy/dt = a
	% note that the constant for the midpoint scheme is identical to that
	% of the trapezoidal scheme when zmid is approximated as (z0+z1)/2
	C = 0.12;
	emax   = C*dzmax;
	% safety margin is factored in later
	dt_opt = dt*sqrt(tol/emax);

	stat.rmse   = NaN;
	stat.maxe  = emax;
	stat.idmax = idmax;
	stat.dt_opt = dt_opt;

	end

	function r = rfun(z)
		dzmid_dt = obj.dz_dt(t+0.5*dt,0.5*(z0+z));
		r = z - z0 - dt*dzmid_dt;
	end

	function A = jfun(z)
		Jt = obj.jacobian(t+0.5*dt,0.5*(z0+z),true);
		%if (obj.aux.surface_flow)
			A  = (obj.aux.I - (0.5*dt)*Jt);
		%else
		%	A = (obj.aux.I_wo_surface_flow - (0.5*dt)*Jt);
		%end
	end % jfun

end % step_midpoint

