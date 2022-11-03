	
	[p,s] = rietkerk_parameter();

	[b,w,h] = rietkerk_equilibrium(p)

	dx = 1;
	L = 10;
	x  = (0:dx:L)';

	nx = length(x);
	y = [b*ones(nx,1);w*ones(nx,1);h*ones(nx,1)];
	dy_dt = rietkerk1d_dt(0,x,y,p,s)

	
