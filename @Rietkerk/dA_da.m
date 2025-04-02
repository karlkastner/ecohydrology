% 2024-12-12 22:05:38.049672356 +0100
function dz_da = dz_da(obj,t,z)
	[b,w,h] = obj.extract1(z);
	% nn = prod(obj.nx);
	nn = numel(b);
	kb = obj.pmu.kb;
	w0 = obj.pmu.w0;
	dw_da = h.*(b + kb.*w0)./(b + kb);
	dz_da = [sparse([],[],[],nn,nn),
		diag(sparse(dw_da)),
		diag(sparse(-dw_da))
		];
end

