% Sun 12 Nov 13:37:24 CET 2023
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
function dz_dt = dz_dt_react(obj,t,z)
	if (isa(obj.p.k,'function_handle'))
		k = obj.p.k(t);
	else
		a = obj.p.k;
	end
	dz_dt = z.*(1 - z./obj.p.k) - obj.p.c.*z.^obj.p.p./(1 + z.^obj.p.p);
end

