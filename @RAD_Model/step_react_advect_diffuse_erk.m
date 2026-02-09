% Wed 13 Mar 13:34:35 CET 2024
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
function [z, stat] = step_react_advect_diffuse_erk(obj,t,z,dt)
	dx   = obj.dx;
	a    = [obj.pmu.vy; obj.pmu.vx];
	e    = [obj.pmu.ey; obj.pmu.ex];
	nfun = @(z) obj.dz_dt_react(t,z,dt);

	z = step_react_advect_diffuse_erk(dt,dx,a,e,nfun,reshape(z,[obj.nx,obj.nvar]),obj.opt.isreal);
	% TODO error estimate
	if (nargout()>1)
		stat = struct('rmse',NaN,'dt0',dt,'flag',0);
	end
end

