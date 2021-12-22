% 2021-07-05 16:48:44.611627340 +0200
function [A,ih] = jacobian(obj,t,z)
	p = obj.p;
	[b,h,w] = obj.extract1(z);

	U  = p.gb.*w./(w + p.kb).*b;
        In = p.a.*h.*(b + p.k.*p.w0)./(b+p.k);

	A = [ p.Db*obj.D2 + (p.cb.*U./b - p.db).*obj.I, ((p.kb./(w.*(p.kb + w))).*U.*p.cb).*obj.I,     obj.Z
	      ((p.a.*h)./(b + p.k) - U./b - In./(b + p.k)).*obj.I, p.Dw*obj.D2 - (p.rw + p.kb./(w.*(p.kb + w)).*U).*obj.I, (In./h).*obj.I
	      (In - (p.a.*h))./(b + p.k).*obj.I,    obj.Z, p.eh*obj.D2 + p.Dh*obj.D1 - (In./h).*obj.I
	];
	ih = [ (b.*p.cb*p.gb.*h.*p.kb)./(h + p.kb).^2
	       b.*( -p.gb.*h.*p.kb./(h + p.kb).^2 + p.a.*w.*p.k.*(1 - p.w0)./(b + p.k).^2 )
              -p.R - p.a.*b.*p.k.*w.*(1 - p.w0)./(b + p.k).^2];
end

