% Mon 31 May 20:20:46 CEST 2021
% Karl Kästner, Berlin
%
%% coefficients of the time-derivative of the Rietkerk-pde
%
% function dz_dt = dz_dt(obj,t,z)
% p : parameter vector
% s : standard deviation of paramter
% TODO variation of dH is not well implemented good -> define at interfaces between cells
function c = dz_dt_coefficient(obj,t,z)
	p = obj.p;
	if (size(z,2)>1)
		[b,w,h] = obj.extract2(z);
	else
		[b,w,h] = obj.extract1(z);
	end
	if (~isvector(z))
		b=b';
		w=w';
		h=h';
	end
		
	n = prod(obj.n);
	if (isa(p.R,'function_handle'))
		R = p.R(t);
	else
		R = p.R;
	end
	
	% uptake of water by plants
	U = p.gb.*w./(w + p.kw).*b;

	% infiltration of water into soil
%	I = a.*h.*(b + p.kb.*p.w0)./(b+p.kb);
	In = p.a.*obj.infiltration_enhancement(b).*h;


if (0)
	o = ones(size(b));
	zz = zeros(size(b));

	% ode coefficients
	cb = [p.cb*p.gb.*w./(w + p.kw) - p.db, zz, p.eb.*o,  zz];
	%cw = [-U./w - p.rw, zz, p.ew.*o, (p.a.*h.*(b + p.kb.*p.w0))./(b + p.kb)];
	% note: it is tempting to write U/w or In/h, but this fails where w or h are
	cw = [-p.rw - (b.*p.gb)./(p.kw + w), zz, p.ew.*o, (p.a.*h.*(b + p.kb.*p.w0))./(b + p.kb)];
	ch = [-(p.a.*(b + p.kb.*p.w0))./(b + p.kb), p.vh.*o, p.eh.*o,   R.*o];

	% stack output
	c = [cb;cw;ch];
end

	c = cell(3,4);

	c(1,1:4) = {p.cb*p.gb.*w./(w + p.kw) - p.db, 0, p.eb,  0};
	c(2,1:4) = {-p.rw - (b.*p.gb)./(p.kw + w), 0, p.ew, (p.a.*h.*(b + p.kb.*p.w0))./(b + p.kb)};
	c(3,1:4) = {-(p.a.*(b + p.kb.*p.w0))./(b + p.kb), p.vh, p.eh,   R};


	% twice any if 2d data passed
%	assert(any(fat(c))==0,sprintf('nan at %g',t));
end % dz_dt

