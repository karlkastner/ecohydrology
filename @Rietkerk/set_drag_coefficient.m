% 2025-08-05 10:44:16.497474970 +0200
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
function C=set_drag_coefficient(obj,z)
	nn = prod(obj.nx); 
	if (0 ~= obj.p.kbC)
		C = obj.p.Cb*((obj.p.kbC+z(1:nn))./(obj.p.kbC+(obj.p.Cb/obj.p.Cv)*z(1:nn))); 
		% extend, TODO fix for 2D
		C = [C(end);C;C(1)];
	else
		C = obj.p.Cv;
	end
	obj.aux.zero_inertia.C = C;
end

