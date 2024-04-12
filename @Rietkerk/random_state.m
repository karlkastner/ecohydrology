% Mon 31 May 20:20:46 CEST 2021
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
%% generate random initial state
%
function [z0] = random_state(obj,mode,R0,mb0,sb0)
	if (nargin()<2)
		mode =2;
	end
	if (nargin()<2)
		R0 = obj.pmu.R;
	end
	if (nargin()<3)
		mb0 = 1;
	end
	if (nargin()<4)
		sb0 = 1;
	end

	switch (mode)
	case {1}
		o = ones(obj.nx,1);
		% note, this does not necessarily work, as the vegetated
		% state can be zero or negative, when R is to small
		% unvegetated state
		[b0,w0,h0] = obj.homogeneous_state([],0);
		z0 = [b0.*o;w0.*o;h0.*o];
		% homogeneously vegetated state
		p = obj.pmu;
		% since it is oscillating between bare and
		p.R = 2*p.R;
		[b1,w1,h1] = obj.homogeneous_state(p,1);
		y1 = [b1.*o;w1.*o;h1.*o];
		p = rand(3*obj.nx,1);
		yr = p.*z0 + (1-p).*y1;
%		z0 = [p(:,1).*b0 + (1-p(:,1)).*b1;
%		      p(:,2).*w0 + (1-p(:,2)).*w1;
%		      p(:,3).*h0 + (1-p(:,3)).*h1];
	case {2} % randomly perturb rainfall
		% this is very noise and
		% the resulting pattern in the unperturbed case becomes labyrinthic
		n = obj.nx;
		if (1 == length(n)) n(2) = 1; end
		p   = obj.pmu;
		p.R        = gamrnd(R0*mb0,R0*sb0,prod(n),1);
		% state pendendent on local water availablity
		[b,w,h] = obj.homogeneous_state(p,2);
		z0  = [flat(b);flat(w);flat(h)];
	case {3} % randomly perturb biomass
		p   = obj.pmu;
		% the initial R must be 1 or higher, otherwise the biomass is decaying to zero
		n   = obj.nx;
		if (1 == length(n)) n(2) = 1; end
		% vegetated state
		p.R = R0;
		[b0, w0, h0] = obj.homogeneous_state(p,1);
		% small amount of biomass (square root of machine precission)
		% interesting, sqrt(eps,single) ~ 3.45e-4 is too small 
		b0_mu = mb0*b0;
		% account for spatial resolution
		b0_sd = sb0*b0*sqrt(1/prod(obj.dx));
		[ap,bp] = gamma_moment2par(b0_mu,b0_sd);
		b       = gamrnd(ap,bp,n);
		w = w0.*ones(n);
		h = h0.*ones(n);
		z0 = [flat(b);flat(w);flat(h)];
	case {4} % bimodal distribution, unvegetated
		n   = obj.nx;
		if (1 == length(n)) n(2) = 1; end
		[b0, w0, h0] = obj.homogeneous_state(obj.pmu,0);
		% vegetated state
		b1 = sb0*15;
		p1 = 0.3/sb0;
		b = repmat(b0,n);
		w = repmat(w0,n);
		h = repmat(h0,n);
		id = randi(prod(n),round(p1*prod(n)),1);
		b(id) = b1;
		z0 = [flat(b);flat(w);flat(h)];
	end
end

