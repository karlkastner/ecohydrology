% Mon 31 May 20:20:46 CEST 2021
% Karl Kästner, Berlin
%
%% time-derivative of the Rietkerk-pde
%
% function dz_dt = dz_dt(obj,t,z)
% p : parameter vector
% s : standard deviation of paramter
% TODO variation of dH is not well implemented good -> define at interfaces between cells
function dz_dt = dz_dt(obj,t,z)
	p = obj.p;
	s = obj.pst;
%%	if (size(z,2)>1)
%		[b,w,h] = obj.extract2(z);
%	else
		[b,w,h] = obj.extract1(z);
%	end
	if (isvector(z))
		b = b';
		w = w';
		h = h';
	end
		
	n = prod(obj.n);

	if (isa(obj.p.R,'function_handle'))
		R = obj.p.R(t);
	end
	
	% uptake of water by plants
	U = p.gb.*w./(w + p.kw).*b;

	% infiltration of water into soil
%	I = a.*h.*(b + p.kb.*p.w0)./(b+p.kb);
	In = p.a.*obj.infiltration_enhancement(b).*h;

	db_dt = p.cb.*U - p.db.*b + p.eb(1).*(obj.aux.D2*b);
	% + (obj.D1c*eb).*(obj.D1c*b);

	%db_dt = cb*U - db_.*b + eb_.*(D2*b);
	dw_dt = In - U - p.rw.*w + p.ew(1).*(obj.aux.D2*w);

	% constant velocity and diffusion
	dh_dt = p.R - In + p.vh(1).*(obj.aux.D1x*h) + p.eh(1)*(obj.aux.D2x*h);

	if (obj.ndim > 1)
		dh_dt = dh_dt + p.vh(2).*(obj.aux.D1y*h) + p.eh(2)*(obj.aux.D2y*h);
	end

	% stack output
	dz_dt = [db_dt; dw_dt; dh_dt];

	% twice any for 2d data
	assert(any(any(isnan(dz_dt)))==0,sprintf('nan at %g',t));
end % dz_dt

