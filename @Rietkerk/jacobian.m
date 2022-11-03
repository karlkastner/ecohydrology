% 2021-07-05 16:48:44.611627340 +0200
% Karl Kästner, Berlin
%
%% jacobian of the Rietkerk model
%
function [A,res] = jacobian(obj,t,z)
	p = obj.p;
	[b,w,h] = obj.extract1(z);

	U  = p.gb.*w./(w + p.kw).*b;
        In = p.a.*h.*(b + p.kb.*p.w0)./(b+p.kb);

	A = [ p.eb*obj.D2 + (p.cb.*U./b - p.db).*obj.I, ((b.*p.cb*p.gb*p.kw)./(p.kw + w).^2).*obj.I,     obj.Z
		((p.a.*h)./(b + p.kb) - (p.gb.*w)./(p.kw + w) - (p.a.*h.*(b + p.kb.*p.w0))./(b + p.kb).^2).*obj.I, ...
			p.ew*obj.D2 + (-p.rw  - (b.*p.gb)./(p.kw + w) + (b.*p.gb.*w)/(p.kw + w).^2).*obj.I, ...
			(In./h).*obj.I;
		(p.a.*h.*p.kb.*(p.w0 - 1))./(b + p.kb).^2.*obj.I,    obj.Z, p.eh*obj.D2 + p.vh*obj.D1 - (In./h).*obj.I
	];
	if (nargout()>1)
	res = [ (b.*p.cb*p.gb.*h.*p.kb)./(h + p.kb).^2
	       b.*( -p.gb.*h.*p.kb./(h + p.kb).^2 + p.a.*w.*p.kw.*(1 - p.w0)./(b + p.kw).^2 )
              -p.R - p.a.*b.*p.kw.*w.*(1 - p.w0)./(b + p.kw).^2];
	end
end

