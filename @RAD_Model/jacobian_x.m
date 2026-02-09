% 2025-07-29 09:55:18.707687971 +0200
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
function J = jacobian_x(obj,t,z,p_react,transpose)
	J = ( obj.aux.Ax ...
	    + p_react*obj.jacobian_react(t,z) ...
	    );
	if (~transpose)
		J = J';
	end
end

