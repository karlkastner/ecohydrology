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
function [z,b,w,h,J,v,e] = single_patch(obj,t,p,bmax,sp)
	nx = obj.nx;
	if (1 == obj.ndim)
		nx = [nx,1];
	end
	if (nargin()<2)
		t = 0;
	end
	if (nargin()<3)
		p = obj.pmu;
	end
	if (nargin()<4)
		bmax = 20;
	end
	if (nargin()<5)
		sp = 1;
	end

	p    = obj.pmu;
	p.R  = 1;
	p.db = 0.25;
	p.gb = 0.05;
	[b0,w,h] = obj.homogeneous_state(t,p);

	if (isscalar(w))
		w = w*ones(nx);
	end
	if (isscalar(h))
		h = h*ones(nx);
	end

	L = obj.L;
	%r = hypot(cvec(x)-L(1)/2,rvec(y)-L(2)/2);
	switch (obj.ndim)
	case {1}
		[x] = obj.x;
		b    = b0 + bmax*normpdf(cvec(x),L(1)/2,sp)/normpdf(0,0,sp);
	case {2}
		[x,y] = obj.x;
		b    = b0 + bmax*normpdf(cvec(x),L(1)/2,sp)*normpdf(rvec(y),L(2)/2,sp)/normpdf(0,0,sp)^2;
	end
	z = [flat(b); flat(w); flat(h)];
end

