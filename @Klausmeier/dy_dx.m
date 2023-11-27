% Wed 30 Jun 12:01:49 CEST 2021
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
%% (c - nu)*D^3*b + (nu*c - c^2 - b^2 - 1)*D^2*b + (c + c*b^2 - c*d + d*nu)*D*b + (d*b^2 - r*b + d)*b
function [dy_dx,A] = dy_dx(x,y)
	b = y(1);
	A = [0, 1;
             -obj.r*b/(b^2 + 1) + obj.d, obj.c];
	dy_dx = A*y(:); 
%	dg_dx = c*g - r*b^2/(b^2 + 1) + d*b;
%	dy_dx = [db_dx;
%                 dg_dx];
if (0)
	b = y(1);
	g = y(2);
	db_dx = g;                                                                   
	dg_dx = c*g - r*b^2/(b^2 + 1) + d*b;
	dy_dx = [db_dx;
                 dg_dx];

end
end

