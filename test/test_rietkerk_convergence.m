	figure(2)
	clf
	meta = vegetation_metadata();
	rkmap = Rietkerk_Map('path_str','mat/');
	rkmap.init();
%	vh = [4,8,12,16];
	vh = 3:18;
	m = length(vh)
	for idx=1:m
		param = {  'L', meta.L, ...
				 'n', round(meta.L/meta.dx), ...
				 'T', meta.To, ...
				 'nt', round(meta.To/meta.dt),  ...
				 ... 'pmu.R', meta.pmu.R, ...
				 'initial_condition', meta.initial_condition, ...
			         'pss.a', 0, ...
				 'pmu.eh' , 1, ...
				 'pmu.vh' , vh(idx)};
		[t,y,rk] = rkmap.run(param{:});
	b = rk.extract1(y);

	% convergence
	S = periodogram(b'-0*mean(b'),rk.L);
	fx = fourier_axis(rk.L,rk.n);
	fdx=fx>0;
	fmu=wmean(S(fdx,end),fx(fdx));
	lmu = fmu*rk.L;
	% n.b.: for 1./1 this is gradual
	nb = round(1/1*rk.L/lmu)

	Sb = periodogram_bartlett(b'-mean(b'),rk.L,nb,rk.n);
	%S_=trifilt1(trifilt1(S,21),21);
	nt=length(t);
	rmse = rms(S-S(:,end))./rms(S(:,end)).*(nt-(1:nt))./(nt-(1:nt)-1);

	subplot(3,m,idx)
	semilogy(t,rmse);
	[mval,mdx]=max(Sb);
	fmax = abs(fx(mdx));	
	
	subplot(3,m,m+idx)
	plot(t,1./fmax)
	hline(1./fmu)
	subplot(3,m,2*m+idx)
	plot(fx,Sb(:,end))
	xlim([0,3*fmu])
	
	end

