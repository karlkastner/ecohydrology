% 2026-08-07 09:46:15.593958351 +0200
function [b,bmax] = homogeneous_state(obj,p)
	if (nargin()<2)
		p = obj.p;
	end
	% ef*(1 - b)*b - mu*b = 0
	% b = 1-mu
	b = 1 - p.mu./p.ef;
	% TODO = bc
	bmax = 1;
end

