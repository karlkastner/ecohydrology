% Mon 31 May 20:20:46 CEST 2021
function [y0] = random_state(obj)
		o = ones(obj.n,1);
if (0)
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
else
		p = obj.pmu;
		p.R = p.R*gamrnd(1,1,obj.n,1);
		[b0,w0,h0] = obj.homogeneous_state(p,1);
		y0 = [b0.*o;w0.*o;h0.*o];
		y0 = max(0,y0);
end
end

