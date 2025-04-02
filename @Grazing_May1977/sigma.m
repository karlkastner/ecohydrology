% Sat 16 Nov 17:04:53 CET 2024
function sigma = sigma(obj,t,z)
	if (~isempty(obj.ptsrel))
		sigma = obj.ptsrel.*z;
	else
		sigma = [];
	end
end

