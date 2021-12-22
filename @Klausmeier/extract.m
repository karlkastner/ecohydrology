% 2021-06-30 22:04:16.711272411 +0200
function [b,w] = extract(obj,y)
	n = obj.n;
	b = y(1:n);
	w = y(n+1:end);
end
