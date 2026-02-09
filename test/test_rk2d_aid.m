% Mon 30 Jun 17:17:28 CEST 2025

if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;
lw = 1;

	% model parameters
	param        = struct();
	
	param.opt.nonlinear_flow = 1;

	% advection coefficient
	param.pmu.vx = [0, 0, 0];
	param.pmu.vy = [0, 0, 0];

	% diffusion coefficient
	param.pmu.ex = [0.1, 0.1, (1-param.opt.nonlinear_flow)*100];
	param.pmu.ey = [0.1, 0.1, (1-param.opt.nonlinear_flow)*100];

	% precipitation (continuous)
	param.pmu.R  = 0.75;

	% mean infiltration coefficient
	param.pmu.a = 0.2;

	% standard deviation of infiltration coefficient in the limit dx -> 0
	param.pss.a  = 0; % 0.2*0.001;
	param.psl.a  = 1;
	param.psdist.a  = 'geometric-ornstein-uhlenbeck';

	% boundary condition
	param.boundary_condition    = {[0,0,1,0],[0,0,1,0];
                                       [0,0,1,0],[0,0,1,0]};

	% reload values of intermediate time steps
	param.opt.loadfinal = false;

	param.L = [16,16]*3;	
	param.L = 2*16*[1,1]*3;	
%	param.L = [1024,1024];	
%	param.L = [4, 4];	

	% initial condition
%	param.initial_condition = 'obj.random_state()';
	param.initial_condition = 'obj.ic_single_patch()';

	% spatial discretization
	dx = 1;

	% number of grid cells
	param.nx = param.L/dx;

	% final time
	param.T = 1;
	% time step for writing output files
	param.opt.dto = 0.1;
	%inf;
	%param.opt.rms_delta_zo_rel_max = 0.1;

	% time step
	param.opt.adapt_time_step=1;
	param.opt.dt           = 1e-3;
	param.opt.dt_min       = 1/400;
	param.opt.dt_max       = inf;
	param.opt.outer_abstol = 1e-3;
	param.opt.outer_reltol = 1e-3;
	param.opt.inner_q      = 0.5;
	param.opt.dt_max_scale_up   = sqrt(2);
	param.opt.dt_min_scale_down = 0;

	% parameters for nonlinear flow
	if (param.opt.nonlinear_flow)
	param.opt.nonlinear_flow_discretization = 'optimal';

	param.pmu.Chezy  = 10;
	param.pss.Chezy  = 0;
	param.psl.Chezy  = 0;
	param.psdist.Chezy  = 'normal';

	param.pmu.lcd = 0.001;
	param.pss.lcd = 0;
	param.psl.lcd = 0;
	param.psdist.lcd   = 'normal';

	param.pmu.zb = zeros(param.nx+2);
	param.pss.zb = 0; 
	param.psl.zb = 0;
	param.psdist.zb = [];
	end

	% data type of output file
	param.opt.compute_class = @double;
	%param.opt.output_class  = @double;
	param.opt.output_class  = @single;
	% output directory
	param.opt.path_str = 'mat/aid/';
	% solver
	%param.opt.inner_solver  = 'gauss-newton';
	param.opt.linear_solver = 'bicgstabl';
	param.opt.inner_q = 0.5;
	param.opt.inner_maxiter = 200;
	param.opt.preconditioner = 'ilu';
	% m to mm
	param.opt.rng = 0;
		param.opt.zero_inertia.input_factor  = 1e-3;
		param.opt.zero_inertia.output_factor = 0.01*86400;
		param.opt.zero_inertia.analytic_jacobian = 0;
	
	

%	if (size(y,2)~=2)
%		y = y';
%	end
%	tY = 365;

	solver       = {'solve_implicit','solve_split','step_aid'}
	% 'solve_split',
	inner_solver = {'gauss-newton','step_advect_diffuse_aid',[]};
 	% 'step_advect_diffuse_implicit_q_fft'
	lsolver = {'bicgstabl','','mldivide'};
	lsolver = {'bicgstabl','','ilu-ks'};
	% ''
	%'ilu-ks'}
	bb = [];

	%id = [1,3];
	id = [1:3];
	solver = solver(id);
	inner_solver = inner_solver(id);
	lsolver = lsolver(id);

	n_step = [];
	runtime = [];
	for idx=1:length(solver)
		param.opt.solver        = solver{idx};
		param.opt.inner_solver  = inner_solver{idx}; 
		param.opt.linear_solver = lsolver{idx};
		rad = Rietkerk(param);
		rad.hashfield_C{end+1} = 'opt.nonlinear_flow';
		rad.hashfield_C{end+1} = 'opt.linear_solver';
		rad.hashfield_C{end+1} = 'opt.inner_solver';
		rad.hashfield_C{end+1} = 'opt.inner_q';
		rad.hashfield_C{end+1} = 'opt.outer_abstol';
		rad.hashfield_C{end+1} = 'opt.outer_reltol';
		rad.hashfield_C{end+1} = 'opt.dt_max_scale_up';
		rad.hashfield_C{end+1} = 'opt.dt_max';
		if (param.opt.nonlinear_flow)
		%rad.hashfield_C{end+1} = 'opt.analytic_jacobian';
		end

		[t,y,out]	  = rad.run();
		b = rad.extract2(y(end,:));
		subplot(2,3,idx)
		imagesc(b)
		bb(:,:,idx) = b;
		out
		if (1 == idx)
			t1 = t;
			y1 = y;
		else
			y1i = interp1(t1,y1,t,'linear');
		end
		n_step(idx) = out.n_step(end);
		runtime(idx) = out.runtime(end)
		%$b = rad.extract2(y(end,:)); imagesc(b); b_fft = b;
	end
	d = bb - bb(:,:,1);
	rmsd = flat(rms(rms(d,2),1))
if (0)
	subplot(2,3,4)
	imagesc(bb(:,:,1)-bb(:,:,2))
	axis equal;
	axis tight
	colorbar
	dydt = rms(diff(y),2)./diff(cvec(t));

	subplot(2,3,5);
	cla
 	semilogy(mid(t),dydt)
	hold on
	plot(t,rms(y-y(end,:),2));

	d = rms(y1i-y,2);
	subplot(2,3,5)
	plot(t,d);
end

subplot(2,3,5); imagesc((bb(:,:,2)-bb(:,:,1))/max(bb(:,:,1),[],'all')); colorbar;
