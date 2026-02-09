
	param     = struct();
	dx        = 1;
	param.L   = 1024;
	param.L   = 128;
	param.nx  = param.L/dx;
	param.opt.linear_infiltration = true;1	
	param.opt.solver   = 'solve_implicit';
	param.opt.inner_solver = 'gauss-newton';
	param.opt.preconditioner = 'ilu';

param.pmu.kgb = 10;
param.pmu.fgb = 0.04;
param.pmu.bevap = 5;
param.pmu.revap = 0.1;
param.pmu.wsat = 50;
param.pmu.rsat = 0.1;

	rad = Rietkerk(param);
	rad.init();

	z   = 20*rand(3*param.nx,1);
	J   = rad.jacobian_react(0,z);

	J_num = jacobian_numeric(@(z) rad.dz_dt_react(0,z),z);

	imagesc((J-J_num)/rms(flat(J)));
	colorbar

	rms(flat(J-J_num))/rms(flat(J))
