% 2024-12-12 22:05:38.049672356 +0100
function d2A_dadz = d2A_dadz(obj,t,z,da)
	[b,w,h] = obj.extract1(z);
	nn = prod(obj.nx);
	kb = obj.pmu.kb;
	w0 = obj.pmu.w0;

	d2w_dadb = h./(b + kb) - (h.*(b + kb*w0))./(b + kb).^2;
	d2w_dadh = (b + kb.*w0)./(b + kb);

	d2A_dadz = [sparse([],[],[],nn,3*nn),
		 diag(sparse(d2w_dadb.*da)),sparse([],[],[],nn,nn),diag(sparse(d2w_dadh.*da))
		 diag(sparse(-d2w_dadb.*da)),sparse([],[],[],nn,nn),diag(sparse(-d2w_dadh.*da))
		];
end

