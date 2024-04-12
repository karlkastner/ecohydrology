% Sun 12 Mar 10:17:12 CET 2023
% nb. this is for 1D models
function [b,w,h] = initial_condition_periodic(obj,fc,a,s,c)
	if (0)
	% the simplest way for a consistent condition is to let the foricing
	% term R oscillate and ignore the diffusion and advection terms
	p = obj.pmu;
	p.R = (0.5*p.R).*(1 + cos(2*pi*fc*obj.x));
	
	[b,w,h] = obj.homogeneous_state(p);
	else
		% TODO 2d
		z = [];
		for idx=1:obj.nvar
			z(:,idx) = a(idx) + s(idx)*sin(2*pi*fc*obj.x) + c(idx)*cos(2*pi*fc*obj.x);
		end
		z = z(:);
	end
end

