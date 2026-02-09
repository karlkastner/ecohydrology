% 2025-04-16 18:52:12.877620704 +0200

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

	%type = 'iso';
	type = 'aniso';

	switch (type)
	case {'iso'}
		L = 512;
		S = 0;
		dx = 1;
		param.opt.dto = 1e2*365;
		param.opt.rms_delta_zo_rel_max = 0.1;
		param.T  = 1e3*365;
	case {'aniso'}
		dx = 2;
		L  = 1024;
		% maximum slope
		S = 1e-3;
		param.opt.dto = 365;
		param.opt.rms_delta_zo_rel_max = inf;
		%param.opt.rms_delta_zo_rel_max = 0.2;
		param.T  = 1e2*365;
	end
	% domain size
	param.L = L;
	nx = param.L/dx;
	n=nx;
	lc = 70;

	x = ((0:nx-1)'+0.5)*dx;
	xe = ((-1:nx)'+0.5)*dx;

	% constant bed slope at outflow boundary = uniform inflow
	dzb_dx = S*cvec(tukeywin(nx+2));
	zb = cumsum(dzb_dx).*dx;

	% mm to m
	hscale = 1e-3;
	param.opt.input_factor = hscale;
	% stretch events over dry periods
	% seconds to days
	param.opt.output_factor = 0.01*86400;

	% zero bed slope at the end = no flow across boundary with neumann condition
if (1)
	zb(end) = zb(end-1);
	zb(1) = zb(2);
end
	s = 1;
		param.opt.inner_maxiter = 80;

	% precipitation
	param.pmu.R  = 0.75;

	% mean infiltration coefficient
	param.pmu.a = 0.2;

	% spatial heterogeneity of parameter a
	param.pss.a  = 0;
	param.psdist.a  = 'lognormal';

	% boundary condition
	param.boundary_condition = {[0,1,0,0],[0,1,0,0]};

	% reload values of intermediate time steps
	param.opt.loadfinal = false;

	% number of grid cells
	param.nx = param.L/dx;

	% final time
	%param.T = 1;	

	% time step
	param.opt.adapt_time_step=1;
	param.opt.dt           = 1/400;
	param.opt.dt_min       = 1/400;
	param.opt.dt_max       = inf;
	param.opt.outer_abstol = 1e-3;
	param.opt.outer_reltol = 1e-3;
	param.opt.inner_q      = 0.5;
	param.opt.dt_max_scale_up   = 1.2; %sqrt(2);
	param.opt.dt_min_scale_down = 0;
	param.opt.analytic_jacobian = true;

	% parameters for nonlinear flow
	if (param.opt.nonlinear_flow)
	param.opt.nonlinear_flow_discretization = 'optimal';

	param.pmu.Chezy  = 5;
	param.pss.Chezy  = 0;
	param.psl.Chezy  = 0; 
	param.psdist.Chezy  = 'normal'; 

	param.pmu.lcd = 0.01;
	param.pss.lcd = 0; 
	param.psl.lcd = 0; 
	param.psdist.lcd = 'normal'; 

	param.pmu.zb = zeros(nx+2,1);
	param.pss.zb = 0; 
	param.psl.zb = 0;
	param.psdist.zb = [];
	end

	% data type of output file
	param.opt.compute_class = @double;
	param.opt.output_class = @double;
	% output directory
	param.opt.path_str = 'mat/nonlinear-flow/';

	% solver
	param.opt.solver   = 'solve_implicit';
	param.opt.inner_solver = 'gauss-newton';
	param.opt.linear_solver = 'mldivide';
	param.opt.input_factor = 1;
	param.opt.output_factor = 1;



	% initial condition
	param.initial_condition = 'obj.random_state()';
	%param.initial_condition = 'obj.ic_single_patch()';
	%param.initial_condition = flat([ 10*(1-cos(2*pi*x/lc)),3*ones(nx,1),14*ones(nx,1)]);
	%rng(0);
	%param.initial_condition = flat([ 10*rand(nx,1),3*ones(nx,1),14*ones(nx,1)]);


	rad = Rietkerk(param);
	rad.p.zb = zb;
	[t,y,out]	  = rad.run();
	if (size(y,2)~=2)
		y = y';
	end

 rad.init_solve();
 h_ = y((2*n+1):(3*n),end);
 [b,w,h__] = rad.extract1(y(:,end));
% [h,u,b,S,a,d] = rad.aux.zero_inertia.interface_values(h__);


y = real(y);
[bb,ww,hh] = rad.extract2(y');
hh = hh';
bb = bb';
ww = ww';
hi = [];ui=[];ai=[];di=[];
for idx=1:length(t)
	[qi,hi(:,idx),he(:,idx),ui(:,idx),bi,Si,ai(:,idx),di(:,idx)] = rad.aux.zero_inertia.interface_values(hh(:,idx),true);
end

figure(1);
clf
%subplot(2,3,1)
%imagesc(y)
subplot(2,3,1)
imagesc(t/365,x,bb)
xlabel('t/y')
subplot(2,3,2)
imagesc(hi)
colorbar
subplot(2,3,3)
imagesc(ui)
colorbar
subplot(2,3,4)
imagesc(ai)
colorbar
subplot(2,3,5)
imagesc(di)
colorbar

%figure(2);
%clf
%subplot(2,3,1)
splitfigure([2,3],[2,1],fflag);
%imagesc(y)
%imagesc(bb)
plot(x,bb(:,end))
ylim([0,30])
xi = inner2outer(x);
% subplot(2,3,2);
splitfigure([2,3],[2,2],fflag);
 plot(xi,hi(:,end));
 ylabel('h')
 yyaxis right;
 plot(xi,-ui(:,end));
 ylabel('-u')
% subplot(2,3,3);
splitfigure([2,3],[2,3],fflag);
 plot(xi,-ai(:,end),'linewidth',lw);
ylabel('Advection $m$/day','interpreter','latex');
ylim([0,20]);
 yyaxis right;
 plot(xi(2:end-1),di(2:end-1,end),'linewidth',lw)
ylabel('Diffusion $m^2$/day','interpreter','latex');
axis square
xlabel('Position $x$/m','interpreter','latex')
% subplot(2,3,4)
splitfigure([2,3],[2,4],fflag);
 plot(xe,zb) 

figure(3)
%clf
subplot(2,3,1)
plot(t,mean(bb))
xlabel('t');
hold on
subplot(2,3,2)
plot(t,mean(hi))
xlabel('t');
hold on
subplot(2,3,3)
plot(t,rms(ui))
subplot(2,2,3)
plot(t,rms(ai))
xlabel('t');
subplot(2,2,4)
plot(t,nanrms(di))
xlabel('t');
hold on

%figure(4);
%clf;
%subplot(2,2,1)
splitfigure([2,2],[4,1],fflag);
cla
loglog(t(1:end-1)/365,out.dt(1:end-1),'.-');
xlabel('Time $t$/year','interpreter','latex');
ylabel('Time step $\Delta t$/day','interpreter','latex');
axis square

%subplot(2,2,2)
splitfigure([2,2],[4,2],fflag);
cla
plot(x,flipud([ ...
100*zb(2:end-1) ...
+ ww(:,end)*[0,1,0,0] ...
+ hh(:,end)*[0,0,1,0] ...
+ bb(:,end).*[1,0,0,0]]),'linewidth',1);
set(gca,'colororder', ...
[0,0.5,0;
0.5,0.25,0.5;
0,0,0.8;
0,0,0]);
set(gca,'xlim',[0,L]);
legend('Biomass $b$/(g/m$^2$)','Soil Water $w$/mm','Surface Water $h$/mm','Bed level $z_b$/cm','location','southwest','interpreter','latex');
xlabel('Position x/m')
axis square

ie = rad.infiltration_enhancement(bb(:,end));fdx = (x>0.25*L & x<0.75*L); 
infiltration_ratio_in_mid_section = mean(ie(fdx)*rad.pmu.a.*hh(fdx,end)) / param.pmu.R

if (pflag)
	ps = 3.5;
%	pdfprint(41,['img/nonlinear-flow-',type,'-time-step'],ps);
%	pdfprint(42,['img/nonlinear-flow-',type,'-pattern'],ps);
	pdfprint(23,['img/nonlinear-flow-',type,'-advection-diffusion'],ps);
end

