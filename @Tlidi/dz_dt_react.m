% Thu 11 Jun 19:03:05 CEST 2026
function dz_dt = dz_dt_react(obj,t,b)
	isvec = 0;
	% convolve kernels
	switch (obj.ndim)
	case {1}
		bf = ifft(obj.aux.Tf.*fft(obj.p.xf.*b));
		bc = ifft(obj.aux.Tc.*fft(obj.p.xc.*b));
	%	bs = ifft(aux.Ts.*fft(b));
	case {2}
		if (isvector(b))
			isvec= true;
			b = reshape(b,obj.nx);
		end
		bf = ifft2(obj.aux.Tf.*fft2(obj.p.xf.*b));
		bc = ifft2(obj.aux.Tc.*fft2(obj.p.xc.*b));
	%	bs = ifft2(aux.Ts.*fft2(b));
	end

	% note that droguette mora has D*Delta b instead of bs
	% note 1/pi was removed to make consistent with DM
	dz_dt = (  b.*(obj.p.bc-b).*obj.p.ef.*exp(bf) ...
		 - obj.p.mu.*b.*exp(bc) ...
		 ... %+ 0*s*D*(bs - aux.int_Rs*b) ...
		);
	dz_dt = real(dz_dt);
	if (isvec)
		dz_dt = dz_dt(:);
	end
end

