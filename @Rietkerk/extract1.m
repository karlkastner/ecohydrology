% Mon  5 Jul 17:45:42 CEST 2021
% Karl Kästner, Berlin
%
%% extract biomass, soil water and surface water from the combined vector
%
function [b,w,h] = extract1(obj,z)
	n = prod(obj.n);
	if (isvector(z))
		b = z(1:n);	
		w = z(n+1:2*n);
		h = z(2*n+1:end);
	else
		b = z(:,1:n);	
		w = z(:,n+1:2*n);
		h = z(:,2*n+1:end);
	end
end

