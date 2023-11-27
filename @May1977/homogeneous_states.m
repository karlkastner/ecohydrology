% Sun 12 Nov 13:41:52 CET 2023
% TODO add trivial solution 0
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
function [S0,lambda] = homogeneous_state(obj,p);
	c = obj.pmu.c;
	k = obj.pmu.k;
	p = obj.pmu.p;
	% p cannot vary
        l = max([length(c),length(k)]);
	m            = p+2;
	z            = zeros(l,m);
	if (issym(k) || issym(c))
		z=sym(z)
	end

	if (~issym(k))
		k = cvec(k)
	end
	if (~issym(c))
		c = cvec(c)
	end
	if (~issym(p))
		p = cvec(p)
	end

	z(:,1)   =  k;
	z(:,2)   = -1;
	z(:,m-2) = z(:,m-2) - c.*k;
	z(:,m-1) = z(:,m-1) + k;
	z(:,m)   = z(:,m) - 1;
	z = fliplr(z);
	if (~issym(z))
       S0 = zeros(l,m-1);
       for idx=1:l
	       S0(idx,:) = roots(z(idx,:));
       end
	else
		eq = z;
		syms z;
		eq = sum(eq.*z.^fliplr(0:length(eq)-1))
		S0 = solve(eq,z,'maxdegree',3);
	end
       S0 = S0.';
%       lambda = -(b - S0.^(p - 1).*p.*(a + r) + S0.^p.*b.*(p + 1));
       % test
%if (0)
 %      rms(-b*S.^(p+1) + (a+r)*S.^p - b*S + a)
%end
end
	

