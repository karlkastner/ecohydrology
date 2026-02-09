% 2025-10-11 16:45:09.104355811 +0200
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
function [z1,stat] = step_trbdf2(obj,t,z0,dt)
		% trapezoidal step from t to t+g*dt
		g  = 2-sqrt(2);

		dz_dt0 = obj.dz_dt(t,z0);
		resg = z0 + 0.5*g*dt*dz_dt0;
		% TODO is the true correct?
		obj.aux.rfun = @(zg) zg - 0.5*g*dt*obj.dz_dt(t+g*dt,zg) - resg;
		obj.aux.jfun = @(zg) (obj.aux.I - 0.5*g*dt*obj.jacobian(t+g*dt,zg,true));
		% solve nonlinear optimization problem
		[zg,stat] = obj.gauss_newton(z0);


		% bdf2 step from t to t+dt
		res1    = (1/g*zg - (1-g)^2/g*z0);
		obj.aux.rfun = @(z1) (2-g)*z1 - (1-g)*dt*obj.dz_dt(t+dt,z1) - res1;
		obj.aux.jfun = @(z1) ((2-g)*obj.aux.I - (1-g)*dt*obj.jacobian(t+dt,z1,true));
		% solve nonlinear optimization problem
		[z1,stat2] = obj.gauss_newton(zg);
		stat.flag  = [stat.flag, stat2.flag];
		stat.iter  = [stat.iter; stat2.iter];
		stat.linear_solver    = [stat.linear_solver, stat2.linear_solver];

		if (nargout()>1)
		% TODO get this from the last call of resgfun
		dz_dtg = obj.dz_dt(t+g*dt,zg);
		dz_dt1 = obj.dz_dt(t+dt,z1);
	
		tol = obj.time_integration_tolerance(z1);
		
		dz = (1/g*dz_dt0 - 1/(g*(1-g))*dz_dtg + 1/(1-g)*dz_dt1);
		[emax,idmax] = max(abs(dz));
	
		C = -0.04*(-2*dt);
		emax = C*emax;
		dt_opt = dt*sqrt(tol/emax);
	
		stat.rmse   = NaN;
		stat.maxe   = emax;
		stat.idmax  = idmax;
		stat.dt_opt = dt_opt;
		end

end % step_trbdf2

