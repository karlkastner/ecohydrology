% Thu  7 Dec 16:41:17 CET 2023
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
%% function [f,g,A] = gnfun(obj,t,dt,z1,z0)
function [f,g,A] = gnfun(obj,t,dt,z1,z0)
	% z1 - z0 = dt*q*dz1 + dt*(1-q)*dz0
	% n = obj.nvar*prod(obj.nx);
	%n=1;
	q   = obj.opt.inner_q;
	% residual r is equal to gradient g
	g = rfun(t,dt,z1,z0);
	if (0)
		e = sqrt(eps);
		% frechet derivative
		g = (rfun(t,dt,z1+e*r,z0)-r)/e;
	end
	% objective
	f = 0.5*(g'*g);

	% bring in form required for java solver
	% TODO, transposing is slow, directly precompute
	%if (obj.aux.)
%		g = reshape(g,obj.nx(1)*obj.nx(2),3)';
	%end

	if (nargout()>2)
	% note that in some cases both A and aa are required
	if (obj.aux.compute_A)
		[Jt,Jrt] = obj.jacobian(t+dt,z1,g,true);
	%	Ar = obj.aux.I - q*dt*Jrt;
		A = (obj.aux.I - (q*dt)*Jt);
	else
		% dummy
		A = [];
	end
	if (obj.aux.compute_aa)
		obj.jacobian_react(t+dt,z1,g,2,obj.aux.tflag);
		obj.aux.aa = -q*dt*obj.aux.aa;
	for idx=1:obj.nvar
			obj.aux.aa(idx,idx,:) = 1.0+obj.aux.aa(idx,idx,:);
		end
	end
	end % if nargout>2

function r = rfun(t,dt,z1,z0)
	% TODO only compute z0 when q < 1, z1 when q>0
	% TODO dz0_dt is fixed and does not need to be recomputed
	dz0_dt = obj.dz_dt(t,z0);
	dz1_dt = obj.dz_dt(t+dt,z1);
	% (z1 - z0)/dt = q*A(z1) + (1-q)*A(z0)
	%r = z1 - z0 - ((1-q)*dt)*dz0_dt - dt*(q*dz1_dt + rhs);
	r = z1 - z0 - ((1-q)*dt)*dz0_dt - (q*dt)*dz1_dt;
end % rfun

end % fun
