% Wed 18 Oct 08:42:39 CEST 2023
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
function [z0,J] = derive_homogeneous_state()
	syms z
	obj = Dodorico();
	obj.make_symbolic();
	obj.p.f = 0;
	obj.p.g = 0;
	dz_dt = obj.dz_dt_react(0,z);
	z0    = solve(dz_dt,z,'MaxDegree',3);
	J     = diff(dz_dt,z);
end

