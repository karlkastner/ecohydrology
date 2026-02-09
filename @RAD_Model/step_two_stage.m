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
%% two stage step, with richardson extrapolation,
%% order higher order accurate than the underlying scheme
%%
function [z,stat] = step_two_stage(obj,t,z0,dt)
	% first half step
	[z1] = obj.aux.fstep(t,z0,dt/2);
	% second half step
	[z2] = obj.aux.fstep(t+dt/2,z1,dt/2);
	% coarse full step
	[z2_] = obj.aux.fstep(t,z0,dt);
	% error estimate
	% TODO set p according to scheme
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
	% TODO choose constant and power according to scheme
	C = 0.24;
	dt_opt = C*dt*sqrt(tol/emax);
	
	stat.rmse   = NaN;
	stat.maxe   = emax;
	stat.idmax  = idmax;
	stat.dt_opt = dt_opt;
	end
end % step_two_stage

