% 2024-12-09 10:04:04.735857954 +0100
function sigma = sigma(obj,t,z)
	[b,w,h] = obj.extract1(z);
	if (~isempty(obj.ptsrel))
		sigma = [obj.ptsrel.b.*b;
			 obj.ptsrel.w.*w
			 obj.ptsrel.h.*h
			];
	else
		sigma = [];
	end
	if (~isempty(obj.pts))
		nn = prod(obj.bx);
		sigma(1:nn) = sigma(1:nn) + obj.pts.b;
		sigma(nn+1:2*nn) = sigma(nn+1:2*nn) + obj.pts.w;
		sigma(2*nn+1:end) = sigma(2*nn+1:end) + obj.pts.h;
	end
end

