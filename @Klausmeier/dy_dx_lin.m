% Wed 30 Jun 12:01:49 CEST 2021
function [dy_dx,A] = dy_dx(x,y)
	b0 = obj.b0;
	r = obj.r;
	d = obj.d;
	c = obj.v;
	A     = [0, 1;
	         -2*r*b0/(b0^2 + 1) + d, c ];
	dy_dx = A*y(:);         

	if (0)                              
	e = y(1);
	g = y(2);
	de_dx = g;                     
	dg_dx = c*g - r*2*b0*e/(b0^2 + 1) + d*e;
	dy_dx_ = [de_dx;
                 dg_dx];
	rms(dy_dx-dy_dx_)
	end
end

