% Sun 11 Feb 22:23:08 CET 2024

javaaddpath('/home/pia/phd/src/lib/mathematics/linear-algebra/multigrid/');

% convergence

%innersolver_C = {'step_advection_diffusion_trapezoidal','step_diffuse_spectral'};
innersolver_C = { ... % 'step_diffuse_spectral', ...
		 'step_advection_diffusion_trapezoidal', ...
		 'step_advect_diffuse_spectral_q' ...
		};
%solver_C = {'euler_forward','solve_split','solve_implicit'};
%solver_C = {'solve_implicit'};
adapt_time_step = [0,0];
%dt  = [0.005/2,0.5,0.5]

%T_pre  = 0;
%dt_pre = 0.5; 

T    = 1000;
%dt_  = [4,2,1,1/2,1/4,1/16,1/32];
%dt_  = [2,1];
%dt_  = [4,2,1,1/2,1/4,1/16,1/32];
%dt_ = 1;
%dt = [1/100,1/10,1];
%dt = [0.5]; %,0.5]; %,0.2,0.1,0.01] %,1/10,1/100];
%dt = 1/400;
dt = 0.25;
tol = [1e-3,3e-4,1e-4]
tol = 3e-4;
b = [];
d=[];

	L     = 128*[1, 1];
	dx    = [1, 1];
	nx    = L./dx; 

	a     = 0.2;
	p_noise = 2; % 2

	s_a   = 0;

	% random generated seeds for repeated experiments
	param   = struct();
	param.opt.loadfinal = true;
	%param.psdist.a = 'geometric-ornstein-uhlenbeck';
	param.initial_condition = 'obj.random_state()';  
	
	% TODO unify
	param.opt.output_class =  @double;
	param.opt.compute_class = @double;
	param.boundary_condition = {'circular','circular'};
	param.pmu.R = 0.7;
	param.pss.a = 0;

	param.L = L;
        param.nx = nx;
        param.opt.dto = T/100;
        param.pmu.a  = a;
        param.pss.a  = 0.0; %0.5;
	param.psdist.a = 'geometric-pink';
	param.opt.rng = 0;
	
	param.pmu.ex=[0.1,0.1,0];
	param.pmu.vx = [0,0,10];

	param.opt.inner2_tol = sqrt(eps);
	param.opt.adapt_time_step = 0; %adapt_time_step(idx);
        param.opt.path_str = sprintf('mat/convergence/');
	%param.psl.a = dt(idx);
	param.opt.inner_q = 0.5;

	% prerun for 100 s to get smooth initial condition
	%param.opt.dt = dtpre;
	%param.T      = T_pre;
	param.opt.solver = 'solve_split';
	%param.opt.loadfinal = false;
	%rk = Rietkerk(param);
	%[t,z,out] = rk.run();
	%z0 = z(end,:);
	param.initial_condition = 'obj.random_state()';  

for idx=1:length(innersolver_C)
	param.opt.inner_solver = innersolver_C{idx};
for tdx=1:length(tol) %dt)

        % run model
	param.T  = T(1);
        param.opt.dt     = dt; %(tdx);
	param.opt.dt_min = dt; %(tdx);
	param.opt.outer_reltol = tol(tdx);
	rk = Rietkerk(param);
	rk.hashfield_C{end+1} = 'opt.dt';		
	rk.hashfield_C{end+1} = 'opt.adapt_time_step';		
	rk.hashfield_C{end+1} = 'opt.inner_solver';		
	rk.hashfield_C{end+1} = 'opt.outer_reltol';		
	%if (strcmp(solver_C{idx},'solve_implicit'))
	%	% TODO do this automatically
	%	rk.hashfield_C{end+1} = 'opt.inner_q';		
	%	rk.hashfield_C{end+1} = 'opt.inner2_tol';		
	%end
	[t,z,out] = rk.run();
	b(:,tdx,idx) = z(end,1:end/3);
	%b(:,:,tdx,idx) = rk.extract2(z(end,:));
	%d(tdx,idx+3) = out.rmse(end);
	%d(tdx,idx+3) = out.esum(end);

	nt = length(tol);
	figure(1)
	subplot(2,nt,tdx+nt*(idx-1))
	imagesc(real(reshape(b(:,tdx,idx),rk.nx)))
	colorbar
	axis equal
	axis tight
end
%if (1 == idx)
	%b(:,2:length(dt_),idx) = repmat(b(:,1,idx),1,length(dt_)-1);
%end

%d_ = b(:,:,idx) - b(:,end,idx);
%d(1:length(dt_),idx) = rms(d_,1);
 
end
d = b(:,:,1) - b(:,:,2);
rmsd = rms(d,1)
figure(2);
for idx=1:length(dt)
	subplot(2,nt,idx)
	imagesc(reshape(d(:,idx),rk.nx));
	colorbar
	axis equal
	axis tight
end
%b(:,:,3) = NaN;
if (0)
d_ = b(:,:,2)-b(:,:,1);
d(:,4) = rms(d_,1);
d_ = b(:,:,3)-b(:,:,1);
d(:,5) = rms(d_,1);
d_ = b(:,:,3)-b(:,:,2);
d(:,6) = rms(d_,1);
end

%subplot(2,2,1)
%loglog(dt,d,'.-');
%legend(solver_C{1:3},'e-sp','e-i','sp-i')
%subplot(2,2,2)
%bar(d(end-1,1:3))
%legend(solver_C)


%subplot(2,2,idx)
%imagesc(b(:,:,idx))
