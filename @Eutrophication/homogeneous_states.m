% 2022-10-20 10:49:45.933121621 +0200 eutro_critical_points.m
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
% dS/dt = a - b S + r S^p/(1+S^q) + d^2/dx^ S
%      0 = (a - b S)(1 + S^p) + r S^p
%      0 = -b S^(p+1)  + (a+r) S^p - b S + a
%      0 = -b (S0+c*cos(x/L))^(p+1)  + (a+r) (S0+c*cos(x/L))^p - b (S0 + c*cos(x/L)) + a - c*D/L^2*cos(x/L)
function [S0,lambda] = homogenous_state(obj);
	a = obj.pmu.a;
	b = obj.pmu.b;
	r = obj.pmu.r;
	p = obj.pmu.p;
	q = obj.pmu.q;
	% p cannot vary
       l = max([length(a),length(b),length(r)]);
	m            = max(p+1,q+2);
	z            = zeros(l,m);
	z(:,m)       = a;
	z(:,m-1)     = -b;
	z(:,m-(p+0)) = z(:,m-(p+0)) + cvec(a+r);
	z(:,m-(q+1)) = z(:,m-(q+1)) - b;

       S0 = zeros(l,m-1);
       for idx=1:l
	       S0(idx,:) = roots(z(idx,:));
       end
       S0 = S0.';
       lambda = -(b - S0.^(p - 1).*p.*(a + r) + S0.^p.*b.*(p + 1));
end
	

