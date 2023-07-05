% Sun 12 Mar 10:17:12 CET 2023
% nb. this is for 1D models
function [b,w,h] = initial_condition_periodic(obj,fc)
	% the simplest way for a consistent condition is to let the foricing
	% term R oscillate and ignore the diffusion and advection terms
	p = obj.pmu;
	p.R = (0.5*p.R).*(1 + cos(2*pi*fc*obj.x));
	
	[b,w,h] = obj.homogeneous_state(p);
end

