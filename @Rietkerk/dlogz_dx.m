% 2021-07-02 12:42:44.980292936 +0200
% Karl Kästner, Berlin
%
%% Rietkerk pde transformed to set of odes through assuming wave-equations
% 
% 
% log-transformed version
%
% note : this is unstable, as the linear version

function dy_dx = rietkerk_dx_log(x,y,p)
	b = exp(y(1));
	w = exp(y(2));
	h = exp(y(3));
	dlb_dx = y(4);
	dlw_dx = y(5);

	% uptake of water by plants
	U = p.gb.*w./(w + p.kb).*b;

	% infiltration of water into soil
	I = p.a.*h.*(b + p.k.*p.w0)./(b+p.k);

	c = p.c;

	% dlb_dx	
	dy_dx(1,1) = dlb_dx;
	dy_dx(2,1) = dlw_dx;
	dy_dx(3,1) =   1/(c-p.Dh)*1/h*(p.R - I);
	% d2lb_dx2
	% d^2/dx^2 = b*(db_dx)^2 + b*d^2lb/dx^2
	dy_dx(4,1) =  -1/p.Db*1/b*(p.cb*U - p.db*b - c*b*dlb_dx) - dlb_dx.^2;
	dy_dx(5,1) =  -1/p.Dw*1/w*( I - U - p.rw*w - c*w*dlw_dx) - dlw_dx.^2;
end

