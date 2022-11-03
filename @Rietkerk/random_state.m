% Mon 31 May 20:20:46 CEST 2021
% Karl Kästner, Berlin
%
%% generate random initial state
%
function [y0] = random_state(obj,mode)
	if (nargin()<2)
		mode =2;
	end
	switch (mode)
	case {1}
		o = ones(obj.n,1);
		% note, this does not necessarily work, as the vegetated
		% state can be zero or negative, when R is to small
		% unvegetated state
		[b0,w0,h0] = obj.homogeneous_state([],0);
		y0 = [b0.*o;w0.*o;h0.*o];
		% homogeneously vegetated state
		p = obj.pmu;
		% since it is oscillating between bare and
		p.R = 2*p.R;
		[b1,w1,h1] = obj.homogeneous_state(p,1);
		y1 = [b1.*o;w1.*o;h1.*o];
		p = rand(3*obj.n,1);
		yr = p.*y0 + (1-p).*y1;
%		y0 = [p(:,1).*b0 + (1-p(:,1)).*b1;
%		      p(:,2).*w0 + (1-p(:,2)).*w1;
%		      p(:,3).*h0 + (1-p(:,3)).*h1];
	case {2}
		p   = obj.pmu;
		n = obj.n;
		if (1 == length(n)) n(2) = 1; end
		p.R = p.R*gamrnd(1,1,n);
		[b1,w1,h1] = obj.homogeneous_state(p,1);
		[b0,w0,h0] = obj.homogeneous_state(p,0);
		fdx = b1<0;
		b1(fdx) = b0(fdx);
		w1(fdx) = w0(fdx);
		h1(fdx) = h0(fdx);
		y1  = [flat(b1);flat(w1);flat(h1)];
		y0  = y1;
		%y0  = max(0,y0);
	case {3}
		p   = obj.pmu;
		n   = obj.n;
		if (1 == length(n)) n(2) = 1; end
		% unvegetated state
		[b0,w0,h0] = obj.homogeneous_state(p,0);
		%[b1,w1,h1] = obj.homogeneous_state(p,1);
		%fdx = b1<0;
		%b1(fdx) = b0(fdx);
		%w1(fdx) = w0(fdx);
		%h1(fdx) = h0(fdx);
		%y1  = [b1;w1;h1];
		%y0  = y1;
		dx = obj.dx;
		sd    = 1.0*sqrt(1/dx.^obj.ndim);
		% small amount of biomass
		rng(1);
		b = b0 + 1e-2*gamrnd(1/sd^2,sd^2,n);
		w = w0.*ones(n);
		h = h0.*ones(n);
		y0 = [flat(b);flat(w);flat(h)];
	end	

end

