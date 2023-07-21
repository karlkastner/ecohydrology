% Mon 31 May 20:20:46 CEST 2021
% Karl Kästner, Berlin
%
%% coefficients of the time-derivative of the Rietkerk-pde
%
% function c = dz_dt_coefficient(obj,t,z)
function c = dz_dt_coefficient(obj,t,z)
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
		
	%n = prod(obj.n);
	
	% uptake of water by plants U_ = U/(wb)
	U_ = obj.p.gb./(w + obj.p.kw);

	% infiltration of water into soil In_ = I/h
	In_ = obj.p.a.*obj.infiltration_enhancement(b);

	% coefficients for c w h
	c = [obj.p.cb.*U_.*w - obj.p.db;
	     -U_.*b - obj.p.rw;
	     -In_];
end % dz_dt_coefficient_react

