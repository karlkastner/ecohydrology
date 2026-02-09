% Mon  2 May 14:18:38 CEST 2022
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
function [z, stat] = step_aid(obj,t,z,dt)
	rmsz0 = rms(z);
	dzdt0 = obj.dz_dt(t,z);

	% weighing of the reaction part scheme
	% 0 : explicit, 1: implicit euler, 0.5: trapezoidal
	p = 0.5;

	% half step x-implicit, y-explicit, reaction-trapezoidal
	% explicit part
	ze        = z + (0.5*dt)*obj.dz_dt_y(t,z,(1-p));
	%(0.5*dt)*obj.dz_dt_y(t,z);
	resfun    = @(z) z - ze - (0.5*dt)*obj.dz_dt_x(t+0.5*dt,z,p);
	jxfun     = @(z) (obj.aux.I - (0.5*dt)*obj.jacobian_x(t+0.5*dt,z,p,true));
	
	% TODO for precition add again half the reaction part
if (0)
	obj.aux.jfun = @jfun;
	obj.aux.rfun = @rfun;
	dz0_dt = obj.dz_dt(t,z0);
	res0   = z0 + (0.5*dt)*dz0_dt;
	% solve nonlinear optimization problem
	[z,stat] = obj.gauss_newton(z0);
end

	[z,exitflag,rmsr,rmsg,res,iter(1),out] = gauss_newton(...
					resfun, ...
					ze+obj.dz_dt_x(t+0.5*dt,ze,p), ...
					obj.opt.nlsolver.abstol, ...
					obj.opt.nlsolver.maxiter, ...
					obj.opt.linear_solver.name, ...
					jxfun, ...
					obj.opt.linear_solver.tol, ...
					obj.opt.linear_solver.maxiter ...
					);

	% half step x-explicit, y-implicit, reaction-trapezoidal
	%ze = z + (0.5*dt)*obj.dz_dt_x(t+0.5*dt,z);
	ze        = z + (0.5*dt)*obj.dz_dt_x(t,z,1-p);
	resfun    = @(z) z - ze - (0.5*dt)*obj.dz_dt_y(t+dt,z,p);
	jyfun     = @(z) (obj.aux.I - (0.5*dt)*obj.jacobian_y(t+dt,z,p,true));

	[z,exitflag,rmsr,rmsg,res,iter(2),out] = gauss_newton(...
					resfun, ...
					ze+obj.dz_dt_y(t+0.5*dt,z,p), ...
					obj.opt.nlsolver.abstol, ...
					obj.opt.nlsolver.maxiter, ...
					obj.opt.linear_solver.name, ...
					jyfun, ...
					obj.opt.linear_solver.tol, ...
					obj.opt.linear_solver.maxiter ...
					);

	% account for temporal noise
	if (0) % ~isempty(sigma))
		sigma = obj.sigma(t,z); 
		% TODO allow for correlation in time
		dz = sqrt(dt)*sigma.*randn(obj.nvar*prod(obj.nx),1);
		z = (z + dz);
	end

	if (nargout()>1)
	% error estimate (second order)
	dz = (obj.dz_dt(t,z) - dzdt0);
	[dmax,idmax]    = max(abs(dz));
	emax            = 0.25*dt*dmax;
	% optimal time step
	tol    = obj.time_integration_tolerance(z);
	dt_opt = dt*sqrt(tol/emax);

	% TODO flag
	linear_solver=struct('iter',NaN);
	stat = struct('maxe',emax,'dt_opt',dt_opt,'flag',exitflag, ...
		'idmax',idmax,'iter',iter,'linear_solver',linear_solver,'rmse',NaN);
	end

	% TODO transpose?
end

