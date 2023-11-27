% Sun 12 Nov 13:37:56 CET 2023
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
	c = obj.pmu.c;
	k = obj.pmu.k;
	p = obj.pmu.p;
	J = 1 + (c.*p.*z.^(2*p - 1))./(z.^p + 1).^2 - (2*z)./k - (c.*p.*z.^p)./(z.*(z.^p + 1));
end

