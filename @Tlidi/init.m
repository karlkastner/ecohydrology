% Thu 11 Jun 19:03:05 CEST 2026
function init(obj)
	%obj.out = struct();
	init@RAD_Model(obj);

	timer = tic();

	% init the convolution kernels
	x = fourier_axis(obj.nx(1)./obj.L(1),obj.nx(1));
	dx = obj.dx();
	
	switch (obj.ndim)
	case{1}
		[x]=fourier_axis(obj.nx(1)./obj.L(1),obj.nx(1));
		Rf = dx/(sqrt(2*pi)*obj.p.lf(1))*exp(-((x-obj.p.dirf(1))/obj.p.lf(1)).^2/2);
		if (1)
			Rf = Rf/(sum(Rf));
		end
		Sf = fft(Rf);
		obj.aux.Tf = sqrt(Sf);
		Rc = dx/(sqrt(2*pi)*obj.p.lc(1))*exp(-(((x-obj.p.dirc(1))/obj.p.lc(1)).^2)/2);
		if (1)
			Rc = Rc/(sum(Rc));
		end
		Sc = fft(Rc);
		obj.aux.Tc = sqrt(Sc);
		%Rs = 1/(sqrt(2*pi)*obj.ls(1))*exp(-(x/lsx).^2/2);
		%Ss = fft(Rs);
		%obj.aux.Ts = sqrt(Ss);
		%obj.aux.int_Rs = sum(Rs,'all')*dx(1);
	case{2}
		dirf = obj.p.dirf;
		dirc   = obj.p.dirc;
		lf = obj.p.lf;
		lc = obj.p.lc;

		[x,y,r]=fourier_axis_2d(obj.nx./obj.L,obj.nx);
		y = rvec(y);
		Rf = dx(1)*dx(2)/(2*pi*lf(1)*lf(2))*exp(-(((x-dirf(1))/lf(1)).^2 + ((y-dirf(2))/lf(2)).^2)/2);
		Sf = fft2(Rf);
		obj.aux.Tf = sqrt(Sf);
		Rc = dx(1)*dx(2)/(2*pi*lc(1)*lc(2))*exp(-(((x-dirc(1))/lc(1)).^2 + ((y-dirc(2))/lc(2)).^2)/2);
		Sc = fft2(Rc);
		obj.aux.Tc = sqrt(Sc);
		%Rs = exp(-s*abs(r.^2));
		%Rs = 1/(2*pi*lsx^2)*exp(-abs(r.^2)./(2*lsx.^2));
		%Ss = fft2(Rs);
		%obj.aux.Ts = sqrt(Ss);
		%obj.aux.int_Rs = sum(Rs,'all')*dx(1)*dx(2);

	end
	switch (func2str(obj.opt.compute_class))
	case {'single'}
		obj.aux.Tf = single(obj.aux.Tf);
		obj.aux.Tc = single(obj.aux.Tc);
	end

	obj.out.runtime_init = toc(timer);
end % init

