% Mon 14 Apr 17:04:44 CEST 2025
% generate a pattern with non-linear flow

if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;
lw = 1

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
	param.pss.a  = 0;
	param.psdist.a  = 'lognormal';

	% boundary condition
	param.boundary_condition    = {[0,0,1,0],[0,0,1,0]};

	% reload values of intermediate time steps
	param.opt.loadfinal = false;

	% domain size
	param.L  = 256*[1];

	% spatial discretization
	dx = 1;

	% number of grid cells
	param.nx = param.L/dx;

	% final time
	param.T = 10*360;	

	% time step
	param.opt.adapt_time_step=1;
	param.opt.dt           = 1/400;
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

	param.pmu.Chezy    = 10;
	param.pss.Chezy    = 0;
	param.psl.Chezy    = 0;
	param.psdist.Chezy = 'normal';

	param.pmu.lcd      = 0.001;
	param.pss.lcd      = 0;
	param.psl.lcd      = 0;
	param.psdist.lcd   = 'normal';

	param.pmu.zb = zeros(param.nx+2,1);
	param.pss.zb = 0; 
	param.psl.zb = 0;
	param.psdist.zb = [];
	end

	% time step for writing output files
	param.opt.dto = 30;
	% data type of output file
	param.opt.compute_class = @double;
	param.opt.output_class = @double;
	% output directory
	param.opt.path_str = 'mat/nonlinear-flow/';
	% solver
	param.opt.solver   = 'solve_implicit';
	param.opt.inner_solver = 'gauss-newton';
% TODO
	param.opt.linear_solver = 'mldivide';
	param.opt.inner_q       = 0.5;
	param.opt.inner_maxiter = 200;
	param.opt.preconditioner = 'ilu';
	% m to mm
	param.opt.zero_inertia.input_factor = 1e-3;
%	param.opt.output_factor = 1;
	param.opt.zero_inertia.output_factor = 0.01*86400;

	% initial condition
	param.initial_condition = 'obj.ic_single_patch()';

	rad = Rietkerk(param);

	rad.hashfield_C{end+1} = 'opt.nonlinear_flow';
	rad.hashfield_C{end+1} = 'opt.linear_solver';
	rad.hashfield_C{end+1} = 'opt.inner_q';
	rad.hashfield_C{end+1} = 'opt.outer_abstol';
	rad.hashfield_C{end+1} = 'opt.outer_reltol';
	rad.hashfield_C{end+1} = 'opt.dt_max_scale_up';
	rad.hashfield_C{end+1} = 'opt.dt_max';
	rad.hashfield_C{end+1} = 'opt.input_factor';
	rad.hashfield_C{end+1} = 'opt.output_factor';

	[t,y,out]	  = rad.run();
	if (size(y,2)~=2)
		y = y';
	end

%y0 = cvec(y(:,end));
%dt = 5;
%T = 1;

if (0)
dz_dt0 = rad.dz_dt_react(0,y0);


zero_inertia = SWE_Zero_Inertia_1d();
zero_inertia.returnmat = true;                          
zero_inertia.zb = zeros(rad.nx+2,1);
zero_inertia.C  = 0.7; 
zero_inertia.lcd = 0.01; 
zero_inertia.L  = rad.L;                                
zero_inertia.n  = rad.nx; 
%zero_inertia.dh_dt_source = dz_dt0(2*end/3+1:end);
t = (0:dt:T)';
h0 = y0(2*end/3+1:end);
h  = zero_inertia.solve(h0,t);
%h_ = zero_inertia.solve(h0,t);
zero_inertia.discretization = 'central';
s = 10;
s_ = 1;
dh_dt = [zero_inertia.dh_dt(h0), s_^2*zero_inertia.dh_dt(1/s*h0)];

dh_dt(1,1)/dh_dt(1,2)

J = zero_inertia.jacobian(h0);
J_ = s_*zero_inertia.jacobian(h0/s);
fdx = find(J);

J(fdx(1))/J_(fdx(1))


figure(1)
clf
subplot(2,2,1)
plot(dh_dt(:,1)./dh_dt(:,2))
subplot(2,2,2)
plot([diag(J)./diag(J_)])

%plot(dh_dt)
pause

if (0)
clf
plot([h0,h(:,end)]);
pause
end % if 0
end % if 0


%	[t,y,out]	  = rad.run();
%	if (size(y,2)~=2)
%		y = y';
%	end
%	yy(:,idx+1) = y(:,end);

% TODO, cumsum values are _not_ correct
if (1)
iter = cvec([out.stat.iter]);
citer = cvec([out.stat.citer]);
figure(1)
clf
subplot(2,3,1)
plot(mid(t),[out.stat.iter],'.-')
ylabel('GN-iterations at step')
yyaxis right
plot(mid(t),citer)
ylabel('Total GN-iterations');

subplot(2,3,2)
plot(mid(t(2:end)),diff(citer)./cvec(diff(t(2:end))),'.-')
hold on
plot(mid(t),(citer)./cvec(mid(t)))
ylabel('GN-steps / simulated day')
xlabel('t / day');
legend('At time t','Average until time t');

subplot(2,3,3)
t=cvec(t);
plot(t,[cvec(out.dt),t./cvec(out.n_step)])
ylabel('dt');
xlabel('t')
legend('at time t','average until t')

end

subplot(2,2,3)
plot(y(:,[1,end]))

subplot(2,2,4)
dy_dt = rms(diff(y,[],2)./rvec(diff(t)));
semilogy(mid(t),dy_dt);

figure(1e4);
di=[];
b=[];
rad.init_solve();
 for idx=1:length(t);
 [b(:,idx),w,h]=rad.extract1(y(:,idx));
 [qi,hi,he,ui,bi,Si,ai,di(:,idx)] = rad.aux.zero_inertia.interface_values(h);
 end;
 subplot(2,2,1);
 imagesc(log(di));
 subplot(2,2,2);
 imagesc(b);
 subplot(2,2,3);
 semilogy(t,max(di));
 yyaxis right;
 plot(t(2:end),iter)

