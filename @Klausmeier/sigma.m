% Sat 16 Nov 17:04:53 CET 2024
function sigma = sigma(obj,t,z)
	[b,w] = obj.extract1(z);
	if (~isempty(obj.pstrel))
		sigma = [obj.pstrel.b.*b;
			 obj.pstrel.w.*w];
	else
		sigma = [];
	end
	if (~isempty(obj.pst))
		sigma(1:prod(obj.nx)) = sigma(1:prod(obj.nx)) + obj.pst.b;
		sigma(prod(obj.nx)+1:end) = sigma(prod(obj.nx)+1:end) + obj.pst.w;
	end
end

