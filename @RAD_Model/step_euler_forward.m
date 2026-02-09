% Mon 16 Oct 09:40:37 CEST 2023
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
function [z,stat] = step_euler_forward(obj,t,z,dt)
	nargin_ = nargin(obj.dz_dt);
	if (2 == nargin_)
		dz_dt = obj.dz_dt(t,z);
	else
		dz_dt = obj.dz_dt(t,z,dt);
	end
	z = z + dt*dz_dt;
	if (nargout()>1)
		% TODO
		stat = struct('rmse',NaN,'dt0',NaN,'flag',0);
	end
end

