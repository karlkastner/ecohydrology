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
function [dy_dx,A] = dy_dx(x,y)
	b0 = obj.b0;
	r = obj.r;
	d = obj.d;
	c = obj.v;
	A     = [0, 1;
	         -2*r*b0/(b0^2 + 1) + d, c ];
	dy_dx = A*y(:);         

	if (0)                              
	e = y(1);
	g = y(2);
	de_dx = g;                     
	dg_dx = c*g - r*2*b0*e/(b0^2 + 1) + d*e;
	dy_dx_ = [de_dx;
                 dg_dx];
	rms(dy_dx-dy_dx_)
	end
end

