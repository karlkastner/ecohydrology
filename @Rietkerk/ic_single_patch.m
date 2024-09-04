% Sun  1 Sep 09:51:23 CEST 2024
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
%% homogeneous (not necessarily stable) states of the Rietkerk system
%
function [z,b,w,h,J,v,e] = single_patch(obj,bmax,sp)
	[b0,w,h] = obj.homogeneous_state();
	if (isscalar(w))
		w = w*ones(obj.nx);
	end
	if (isscalar(h))
		h = h*ones(obj.nx);
	end
	if (nargin()<2)
		bmax = 20;
	end
	if (nargin()<3)
		sp = 1;
	end

	L = obj.L;
	[x,y] = obj.x;
	%r = hypot(cvec(x)-L(1)/2,rvec(y)-L(2)/2);
	b    = b0 + bmax*normpdf(cvec(x),L(1)/2,sp)*normpdf(rvec(y),L(2)/2,sp)/normpdf(0,0,sp)^2;
	z = [flat(b); flat(w); flat(h)];
end

