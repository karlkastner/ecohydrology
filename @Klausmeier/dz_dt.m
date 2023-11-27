% 2021-06-30 22:04:16.711272411 +0200
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
%% c.f. Klausmeier 1999
function dz_dt = dz_dt(obj,t,z)

	[b,w] = obj.extract1(z);

	uptake = obj.p.g.*w.*b.*b;

	db_dt = obj.p.c.*uptake     - obj.p.d.*b + obj.p.eb(1)*(obj.aux.D2x*b) + obj.noise.b;
	dw_dt = obj.p.r - obj.p.l.*w - uptake   + obj.p.vw(1).*(obj.aux.D1x*w) + obj.p.ew(1)*(obj.aux.D2x*w) + obj.noise.w;
	if (length(obj.nx)>1)
		db_dt = db_dt + obj.p.eb(2)*(obj.aux.D2y*b);
		dw_dt = dw_dt + obj.p.vw(2)*(obj.aux.D1y*w) + obj.p.ew(2)*(obj.aux.D2y*w);
	end

	dz_dt = [db_dt; dw_dt];
end

