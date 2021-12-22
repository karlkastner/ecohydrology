% Mon  5 Jul 17:45:42 CEST 2021
function [b,w,h] = extract2(obj,z)
	[b,w,h] = obj.extract1(z);
	if (length(obj.n)>1)
		b = reshape(b,obj.n);
		w = reshape(w,obj.n);
		h = reshape(h,obj.n);
	end
end

