% Thu  7 Dec 09:29:07 CET 2023
function [J,Jr] = jacobian(obj,t,z,r,transpose)
	Jr = obj.jacobian_react(t,z,r);
	if (transpose)
		J  = Jr + obj.aux.AD;
	else
		J  = Jr' + obj.aux.AD';
		Jr = Jr';
	end
end

