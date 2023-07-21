% Mon 31 May 20:20:46 CEST 2021
% Karl Kästner, Berlin
%
%% coefficients of the time-derivative of the Rietkerk-pde
%
% function c = dz_dt_coefficient_react_homogeneous(obj,t,z)
function c = dz_dt_coefficient_react_inhomogeneous(obj,t,z)
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

	if (isa(obj.p.R,'function_handle'))
		R = obj.p.R(t,obj.x,obj.y);
	else
		R = flat(obj.p.R).*ones(prod(obj.n),1);
	end

	% infiltration of water into soil
	In = obj.p.a.*obj.infiltration_enhancement(b).*h;

	c = [zeros(prod(obj.n),1);
             In;
	     R;
	     ];
end

