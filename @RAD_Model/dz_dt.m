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
function dz_dt = dz_dt(obj,t,z)
	% a = obj.afun(t);
	% dS/dt = a - b S + r S^p/(1+S^p) + e d^2/dx^S
	dz_dt = obj.dz_dt_react(t,z) + obj.aux.AD*z;
	%dS_dt = a - b*S + r*S.^p./(1 + S.^q) + e*reshape(laplacian*S(:),n(1),n(2));
end

