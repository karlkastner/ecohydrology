% 2023-10-16 10:36:20.158690521 +0200
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
function [dz_dt] = dz_dt_react(obj, t, z)
	[b, w] = obj.extract1(z);
	if (isa(obj.p.r,'function_handle'))
		r = obj.p.r(t);
	else
		r = obj.p.r;
		% *(1 + 0.1*randn(obj.nx,1));
	end

	uptake = obj.p.g.*w.*b.*b;
	db_dt = obj.p.c.*uptake - obj.p.d.*b;
	% obj.noise.b;
	% db_dt = db_dt + 1e-1*b.*randn(size(b));
	dw_dt = r - obj.p.l.*w - uptake;
	% + obj.noise.w;
	dz_dt = [db_dt; dw_dt];

end



