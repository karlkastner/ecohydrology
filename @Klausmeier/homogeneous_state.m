% 2021-06-30 13:06:19.787743754 +0200
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
function [b0,w0,J,e,v] = homogeneous_state(obj,p)
%	r = obj.r;
%	d = obj.d;
%	b0 = [
%	                               0,                                              
%	 (r + (r^2 - 4*d^2)^(1/2))/(2*d),                                                
%	 (r - (r^2 - 4*d^2)^(1/2))/(2*d) ];
	if (nargin()<2)
		p = obj.pmu;
	end
	g = p.g;
	l = p.l;
	d = p.d;
	c = p.c;
	r = p.r;

	b0 = [ ((-g*(- g*c^2*r^2 + 4*l*d^2))^(1/2) + c*g*r)/(2*d*g)
	      -((-g*(- g*c^2*r^2 + 4*l*d^2))^(1/2) - c*g*r)/(2*d*g)];
	w0 = d./(b0*c*g);

	for idx=1:2
	J(:,:,idx) = jacobian(b0(idx),w0(idx));
	[v(:,:,1),e_] = eig(J(:,:,idx));
	e(:,idx) =diag(e_);
	end
function J = jacobian(b,w)
	J = [2*b*c*g*w - d,     b^2*c*g; 
                  -2*b*g*w, - g*b^2 - l];
end
end

