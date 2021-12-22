% Wed 30 Jun 12:01:49 CEST 2021
% (c - nu)*D^3*b + (nu*c - c^2 - b^2 - 1)*D^2*b + (c + c*b^2 - c*d + d*nu)*D*b + (d*b^2 - r*b + d)*b
function [dy_dx,A] = dy_dx(x,y)
	b = y(1);
	A = [0, 1;
             -obj.r*b/(b^2 + 1) + obj.d, obj.c];
	dy_dx = A*y(:); 
%	dg_dx = c*g - r*b^2/(b^2 + 1) + d*b;
%	dy_dx = [db_dx;
%                 dg_dx];
if (0)
	b = y(1);
	g = y(2);
	db_dx = g;                                                                   
	dg_dx = c*g - r*b^2/(b^2 + 1) + d*b;
	dy_dx = [db_dx;
                 dg_dx];

end
end

