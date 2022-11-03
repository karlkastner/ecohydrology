% Mon  5 Jul 17:45:42 CEST 2021
% Karl Kästner, Berlin
%
%
%
function [b,w,h] = extract2(obj,z)
	[b,w,h] = obj.extract1(z);
	nt = size(z,1);
	if (length(obj.n)>1)
		b = reshape(b,[nt,obj.n(1),obj.n(2)]);
		w = reshape(w,[nt,obj.n(1),obj.n(2)]);
		h = reshape(h,[nt,obj.n(1),obj.n(2)]);
	end
end

