% Fri 23 May 12:41:04 CEST 2025
% generate a pattern with non-linear flow

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
	%param.pmu.R  = 0.75;
	%of = 0.01;
	T_interval = 20;
	T_duration = 0.05;
	%T_duration = 0.1;
	scale = T_interval/T_duration;
	of = scale*0.01;	
%	Rbar = 0.75;
%	Rbar = 1;
%	param.pmu.R  = @(t,~) scale*Rbar*(mod(t,10)<0.1);
	
	%20*360;
	%T= 360;	
	T=50;
	%T = 20;

	if (T_duration < inf)
	R0 = 1.2;
	R0 = 1.0;
	[t_event,r_event] = generate_rainfall(T,T_interval,T_duration,R0*T_interval/T_duration);
	param.opt.r_event = [0;flat(r_event)];

	% precipitation
	param.pmu.R  = @(t,edx) param.opt.r_event(edx);
	param.opt.tevent = t_event;
	else
	param.pmu.R  = R0;
	param.opt.r_event =[];
	param.opt.t_event = [];
	end

	% mean infiltration coefficient
	param.pmu.a = 0.2*scale;


	% standard deviation of infiltration coefficient in the limit dx -> 0
	param.pss.a  = 0.2*0.001;
	param.psl.a  = 1;
	param.psdist.a  = 'geometric-ornstein-uhlenbeck';

	% boundary condition
	param.boundary_condition    = {[0,0,1,0],[0,0,1,0];
                                       [0,0,1,0],[0,0,1,0]};

	% reload values of intermediate time steps
	param.opt.loadfinal = false;

	% domain size
	% 63
				% T=1	2	4	8			20
				% rt in sec
	param.L  = 63*[1,1];	%   21
%	param.L  = 129*[1,1];	%  162	245	337	706p
%	param.L  = 189*[1,1];	%  351
%	param.L  = 255*[1,1];	%  597		1194p	2827s,136MB (dto 30)	5654p	
%	param.L  = 513*[1,1];	% 2406p		4776p	11308p (dto 360)	23378 (dto 180)
%	param.L  = 1023*[1,1];  %						93512 (p)
%	random ic								20
%	255									2075s,221MB (0.05%)	
%	513									9712,477MB (0.1%)

	% dtmin = 1/400	
	% T=30, Tr=0.1,Ti=10		impl   aid
	% 63			24	21
	% 129			105     68	
	% 255   T10		96	103
	

	% initial condition
	param.initial_condition = 'obj.random_state()';
	param.initial_condition = 'obj.ic_single_patch()';

	% spatial discretization
	dx = 1;

	% number of grid cells
	param.nx = param.L/dx;

	% final time
	param.T = T;
	% time step for writing output files
	param.opt.dto = inf;
	param.opt.rms_delta_zo_rel_max = 0.5;

	% time step
	param.opt.adapt_time_step=1;
	param.opt.dt           = 1/400;
	param.opt.dt_min       = 1e-4; %1/400;
	param.opt.dt_max       = inf;
	param.opt.outer_abstol = 1e-3;
	param.opt.outer_reltol = 1e-3;
	param.opt.inner_q      = 0.5;
	param.opt.dt_max_scale_up   = sqrt(2);
	param.opt.dt_min_scale_down = 0;

	% parameters for nonlinear flow
	if (param.opt.nonlinear_flow)
	param.opt.zero_inertia.discretization = 'optimal';

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
	param.opt.output_class  = @single;
	% output directory
	param.opt.path_str = 'mat/nonlinear-flow/';
	% solver

	if (0)	
	param.opt.inner_solver = 'step_advect_diffuse_aid';
	param.opt.zero_inertia.integration_scheme  = 'aid';
	param.opt.zero_inertia.linear_solver_name  = 'direct'; %ilu-ks';
	%param.opt.zero_inertia.integration_scheme  = 'trapezoidal';
	%param.opt.zero_inertia.linear_solver_name  = 'ilu-ks';
	else 
	param.opt.solver        = 'solve_implicit';
	param.opt.inner_solver  = 'gauss-newton';
	param.opt.linear_solver = 'bicgstabl';
	param.opt.inner_q = 0.5;
	param.opt.inner_maxiter = 200;
	param.opt.preconditioner = 'ilu';
	param.opt.zero_inertia.integration_scheme  = '';
	end

	% mm to m
	param.opt.zero_inertia.input_factor  = 1e-3;
	% 10 days to 0.1 day in seconds
	param.opt.zero_inertia.output_factor       = of*86400;
	if (1)
	else
	end
	param.opt.zero_inertia.analytic_jacobian = 0;
	
	% maybe this should not go into opt
	param.opt.linear_infiltration = 0;
	% 10 mm/h/w0 = 10 mm/h/0.2 to mm/day
	param.pmu.a = 1200;
	param.pmu.db = 0.1;
	
	rad = Rietkerk(param);

	rad.hashfield_C{end+1} = 'opt.nonlinear_flow';
	rad.hashfield_C{end+1} = 'opt.linear_solver';
	rad.hashfield_C{end+1} = 'opt.inner_q';
	rad.hashfield_C{end+1} = 'opt.outer_abstol';
	rad.hashfield_C{end+1} = 'opt.outer_reltol';
	rad.hashfield_C{end+1} = 'opt.dt_max_scale_up';
	rad.hashfield_C{end+1} = 'opt.dt_max';
	rad.hashfield_C{end+1} = 'opt.dt_min';
	rad.hashfield_C{end+1} = 'opt.zero_inertia.integration_scheme';
	rad.hashfield_C{end+1} = 'opt.zero_inertia.analytic_jacobian';
	rad.hashfield_C{end+1} = 'opt.r_event';
	rad.hashfield_C{end+1} = 'opt.linear_infiltration';
	[t,y,out]	  = rad.run();
	if (size(y,2)~=2)
		y = y';
	end
tY = 365;

figure(1);
clf
b = rad.extract2(y(:,end));
subplot(2,3,1)
imagesc(rad.x,rad.x,b)
colorbar
hold on

subplot(2,3,2)
dy_dt = rms(diff(y,[],2))./diff(rvec(t));
dy_dt(2,:) = max(abs(diff(y,[],2)))./diff(rvec(t));
plot(mid(t)/tY,dy_dt)
ylabel('||dy/dt||')
dy = rms(y-y(:,end));
%hold on;
yyaxis right
plot(t/tY,dy./(t(end)-t));
ylabel('||y-y(T)||/(t-T)')
xlabel('Simulated time / years');

subplot(2,3,3)
plot(t(1:end-1)/tY,out.dt(1:end-1),'.')
ylabel('Time step length / days');
%xlabel('Simulated time / years');
yyaxis right
dt_=(diff(t(1:end-1))./diff(out.n_step(1:end-1)));
plot(mid(t(1:end-1))/tY,dt_)
xlabel('Simulated time / years');

subplot(2,3,4);
bb=[];
for idx=1:length(t);
 b=rad.extract2(y(:,idx));
 bb(:,idx) = b(:,round(end/2));
 end;
 imagesc(t,rad.x,bb)

try
subplot(2,3,5)
cla
cgn_iter=[0,out.stat.citer];
gn_iter = [out.stat.iter];
iter = diff(cgn_iter)./diff(out.n_step(1:end));
iter_ = diff([0,out.stat.citer])./diff(t);
plot(mid(t)/tY,iter);
ylim([0,1.05*max(iter)])
ylabel('GN-steps per time step')
yyaxis right;
plot(mid(t)/tY,iter_)
ylabel('GN-steps per simulated day')
xlabel('Simulated time / year');
catch
end

try
subplot(2,3,6)
cla
iter = arrayfun(@(x) sum(x.linear_solver.iter),out.stat(2:end));
%plot(iter(1:end-1));
plot(mid(t(1:end-1))/tY,iter(1:end-1))
ylabel('Linear solver steps per time step')
yyaxis right;
diter_dt = iter(1:end-1)./out.dt(2:end-1);
plot(mid(t(1:end-1))/tY,diter_dt)
i_ = ([gn_iter*10+(iter)]);
di_dt_= i_./out.dt(2:end);
%diff(rvec(t));            
hold on;
plot(mid(t)/tY,di_dt_)  
hline(mean(diter_dt),'color','r','linestyle','--')
ylabel('Linear solver steps per simulated day')
xlabel('Simulated time / year');
catch
end
%y0 = cvec(y(:,end));
%dt = 5;
%T = 1;

figure(2);
clf
db = rad.extract2(y(:,end)-y(:,end-1));
%figure(2);
imagesc(db);
colorbar


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
end

C = [0.7];
yy = cvec(y0);
for idx=1:length(C)

% continue to run with nonlinear flow for a day
%param.pmu.zb    = 0;
param.pmu.zb = zeros(rad.nx+2,1);
param.pss.zb = 0;
param.psdist.zb = 'normal';
param.psl.zb = [];

param.pmu.lcd = 0.01;
param.pss.lcd = 0;
param.psdist.lcd = 'normal';
param.psl.lcd = [];
% W=1;T=1; L=40; h=10; Q=L*h/T; dh = 1; S = dh/L; normal_flow_roughness(Q,h,W,S)
param.pmu.Chezy    = C(idx);
param.pss.Chezy = 0;
param.psdist.Chezy = 'normal';
param.psl.Chezy = [];

%	param.opt.solver   = 'solve_split';
param.opt.solver   = 'solve_implicit';
param.opt.inner_solver = 'gauss-newton';
param.boundary_condition    = {[0,0,1,0],[0,0,1,0]};
%	param.opt.inner_solver = 'step_advect_diffuse_spectral';

if (param.opt.nonlinear_flow)
param.pmu.ey(3) = 0;
param.pmu.vx(3) = 0;
param.pmu.vy(3) = 0;
end
%T = 100;
T = 365*40;
%T =1;
param.T         = T;
param.opt.dto = 30;
%param.opt.dto = 180;
% keep time step constant
	param.opt.adapt_time_step = 0;
	param.opt.dt = dt;
	param.opt.adapt_time_step=1;
	param.opt.dt           = 1/400;
	param.opt.dt_min       = 1/400;
	param.opt.dt_max       = inf;
	param.opt.outer_abstol = 1e-3;
	param.opt.outer_reltol = 1e-3;
	param.opt.dt_max_scale_up   = sqrt(2);
	param.opt.dt_min_scale_down = 0;

param.initial_condition = max(y0,0);
param.opt.inner_maxiter = 200;
%param.opt.zero_inertia.discretization = 'optimal';
param.opt.zero_inertia.discretization = 'upwind';
%param.opt.inner_q = 'midpoint';
%param.opt.inner_q = 'trapezoidal';
param.opt.inner_q = 0.5;
%param.opt.inner_q = 1;
param.opt.loadfinal = false;
%param.opt.zero_inertia.discretization = 'central';

% TODO zi velocity and discharge


% compare wl
	param.opt.preconditioner = 'ilu';

	rad = Rietkerk(param);
	rad.hashfield_C{end+1} = 'opt.nonlinear_flow';
	rad.hashfield_C{end+1} = 'opt.linear_solver';
	rad.hashfield_C{end+1} = 'opt.inner_q';
	rad.hashfield_C{end+1} = 'opt.outer_abstol';
	rad.hashfield_C{end+1} = 'opt.outer_reltol';
	rad.hashfield_C{end+1} = 'opt.dt_max_scale_up';
	rad.hashfield_C{end+1} = 'opt.dt_max';
	[t,y,out]	  = rad.run();
	if (size(y,2)~=2)
		y = y';
	end
	yy(:,idx+1) = y(:,end);
end

if (0)

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
end

end

% plot time step limitation
[b,w,h] = rad.extract1(y'); 
figure(1e4);
 clf;
subplot(2,2,1)
plot(t(2:end),floor(([out.stat.idmax]-1)/prod(rad.nx))+1,'*-');
 set(gca,'ytick',1:3,'yticklabel',{'b','w','h'});
 ylim([0.5,3.5]);
 yyaxis right;
 plot(t,mean(h,2))
subplot(2,2,2)
plot(t(2:end),[out.stat.limiting_part])
set(gca,'ytick',[0,1],'yticklabel',{'reaction','diffusion'})
ylim([-0.5,1.5])
yyaxis right
plot(t(1:end),[out.dt])
