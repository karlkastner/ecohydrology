% 2023-10-16 17:26:48.921654888 +0200
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

		pmu = struct(... 
			 'c',   1          ... % conversion rate of water to biomass (J in Klausmeier)
			,'d',  0.045       ... % death rate of biomass (M in klausmeier)	
			,'eb', 1*[1,1]     ... % diffusion of biomass (D in klausmeier)
			,'ew', [0,100]     ... % diffusion of water (not in klausmeier)
			,'g',   1          ... % water uptake (R in klausmeier)
			,'l',   1          ... % evaporation (L in clausmeier)
			,'r',   0.077      ... % precipitation (A in Klausmeier)
			,'vw', 182.5*[1,0] ... % water runoff velocity
		);

syms b w r;
syms c d g l r;


km = Klausmeier();
km.nx = 1;
km.p.c = c;
km.p.d = d;
km.p.g = g;
km.p.l = l;
km.p.r = r;

dz_dt = km.dz_dt_react(0,[b;w]);

w_=solve(dz_dt(1),w)
b_ = solve(subs(dz_dt(2),w,w_),b)

km = Klausmeier();
km.p = km.pmu;
km.p.r = 0.45;
[b0,w0] = km.homogeneous_state()
km.dz_dt_react(0,[b0(1),w0(1)])
km.dz_dt_react(0,[b0(2),w0(2)])

J = [diff(dz_dt,b), diff(dz_dt,w)]

[v,e] = eig(J); e=diag(e); simplify(e,'ignoreanalyticconstraints',true), e1=subs(subs(e,w,w_(1)),b,b_(1))

