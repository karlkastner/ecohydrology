	meta = vegetation_metadata();
	%end
	pflag = meta.pflag;
	fflag = pflag;
	
	L = meta.example.L/8;
	dx = meta.example.dx;
	sd_a = 0;
	vh   = meta.example.vh;

	base  = sprintf(    'rietkerk1d-L-%g-dx-%g-Ti-%g-To-%g' ...
				,round(L),dx,meta.mopt.Ti,meta.mopt.To);
	oname_ = base;

	rkmap = Rietkerk_Map('path_str','mat/');
	rkmap.init();

	param = struct();
	param.L = L;
	param.n = round(L/dx);
	param.T = meta.mopt.To;
	param.nt = round(meta.mopt.To/meta.mopt.dt);
	param.pss.a = sd_a;
	param.pmu.eh = meta.mopt.eh;
	param.pmu.vh = vh;
	param.initial_condition = 'random';
	[t,y,rk]   = rkmap.run(param);
	y0 = rk.initial_condition_from_central_frequency(y);

	param.initial_condition = y0;
	% run model
	[t,y,rk] = rkmap.run(param);
	y = double(y);

	param.zero_inertia = true;
	param.pmu.dzb_dx   = 0.02;
global h0
if (1)
	y0 = 0.*y0;
	h0 = normpdf(rk.x,mean(rk.x),10)*20;
	h0 = cvec(h0);
end
	y0 = [0.*h0;0.*h0;h0];
	param.initial_condition = y0;
		
	
	[t,y_,rk] = rkmap.run(param);
