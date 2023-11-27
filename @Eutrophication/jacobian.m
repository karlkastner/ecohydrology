% Mon 23 Oct 17:11:21 CEST 2023
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
function J  = jacobian(obj,t,z)
	p = obj.pmu.p;
	q = obj.pmu.q;
	r = obj.pmu.r;
	b = obj.pmu.b;
	a = obj.pmu.a;
	
	J = (z.^(p - 1).*p.*r)./(z.^q + 1) - b - (z.^p.*z.^(q - 1).*q.*r)./(z.^q + 1).^2;
end

