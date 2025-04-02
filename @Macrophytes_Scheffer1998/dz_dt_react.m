% Fri 15 Nov 19:09:40 CET 2024
function dz_dt = dz_dt_react(t,z)
	[e,v] = obj.extract1(z);
	dz_dt [ obj.p.re.*E*(1 - (e./obj.p.e0).*(obj.p.hv + V)./obj.p.hv);
		obj.p.rv.*V*(1 - v.*(obj.he.^obj.p.p - e.^obj.p.p)./obj.p.he.^obj.p.p ...
	      ];
end % dz_dt_react

