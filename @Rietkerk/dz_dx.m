% Fri  2 Jul 10:15:04 CEST 2021
%
% rietkerk model transformed to set of odes
% time eliminated with travelling wave assumption:
%	t = 1/c*x
%
% note : this system does not seem stable, why
%
function dy_dx = rietkerk_dx(x,y,p)
	b = y(1);
	w = y(2);
	h = y(3);
	db_dx = y(4);
	dw_dx = y(5);

	% uptake of water by plants
	U = p.gb.*w./(w + p.kb).*b;

	% infiltration of water into soil
	I = p.a.*h.*(b + p.k.*p.w0)./(b+p.k);

	c = p.c;
	
	dy_dx(1,1) = db_dx;
	dy_dx(2,1) = dw_dx;
	dy_dx(3,1) =   1/(c-p.Dh)*(p.R - I);
	dy_dx(4,1) =  -1/p.Db*(p.cb*U - p.db*b - c*db_dx);
	dy_dx(5,1) =  -1/p.Dw*( I - U - p.rw*w - c*dw_dx);
end


