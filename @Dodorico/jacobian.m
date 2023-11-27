% Tue 24 Oct 09:52:54 CEST 2023
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
function J = jacobian(obj,t,z,p)
	if (nargin()<4)
		p = obj.pmu;
	end
	a = p.a;
	b = p.b;
	ze = p.ze;
	zmax = p.zmax;
	l0 = p.l0;
	w0 = p.w0;
	J = [(b.*w0)./(l0 + b.*z).^2 - a.*(z - zmax) - a.*(z - ze)];
end

