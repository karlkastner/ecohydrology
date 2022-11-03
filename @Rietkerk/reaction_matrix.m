% 2022-05-03 13:09:24.580982964 +0200
function [A,ih] = reaction_matrix(obj,y)
	[b,w,h] = obj.extract1(y);
	p = obj.p;
	z = zeros(size(b));
	o = ones(size(b));
	A = [0.5*(p.cb*p.gb.*w)./(p.kw + w) - p.db, 0.5*b.*(p.cb*p.gb)./(p.kw + w), z;
	     - 0.5*w*p.gb./(p.kw + w) + 0.5*p.a.*h./(b + p.kb), -p.rw - 0.5*(b*p.gb)/(p.kw + w), 0.5*p.a.*b./(b + p.kb) + (p.a.*p.kb.*p.w0)./(b + p.kb) ;
	     -(p.a.*h)./(b + p.kb), z, - p.a.*p.kb*p.w0./(b + p.kb) ];
	ih = [z;z;p.R.*o];
end

