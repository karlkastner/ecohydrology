function init_solve(obj)
	% convert solver to string for switch-case
	if (isa(obj.opt.solver,'function_handle'))
		solver_str = func2str(obj.opt.solver);
	else
		solver_str = obj.opt.solver;
	end

	%
	% init solvers
	%
	tic();
	obj.out = struct();
	obj.aux.compute_S = false;
	switch (solver_str)
	case {'solve_split'}
		% initialize fourier transform of impulse responses of linear
		% advection-diffusion part
		switch (obj.opt.inner_solver)
		case {'step_advect_diffuse_trapezoidal'}
			obj.init_fourier_matrices();
		case {'step_advect_diffuse_spectral'}
			% nothing to do
		case {'step_advect_diffuse_implicit_q_fft'}
			% nothing to do
		otherwise
			error(sprintf('unimplemented inner solver %s',obj.opt.inner_solver));
		end
		obj.aux.fstep = @obj.step_split;
	case {'step_integrating_factor'}
		obj.aux.fstep = @obj.step_integrating_factor;
	case {'step_react_advect_diffuse_erk'}
		obj.aux.fstep = @obj.step_react_advect_diffuse_erk;
	case {'euler_forward'}
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_euler_forward;
	case {'solve_implicit'}
		% TODO this should be avoided for MG
		obj.init_advection_diffusion_matrix();
		switch (obj.opt.innersolver2)
		case {'gmres','bicgstabl'}
			switch (obj.opt.preconditioner)
			case {'multigrid-java'}
				obj.aux.compute_aa = true;
				init_multigrid_java();
				obj.aux.mg_j.nmaxiter = 1;	
			case {'ilu','',[]}
				obj.aux.compute_aa = false;
			end
			obj.aux.compute_A = true;				
		case {'multigrid'}
			obj.mg_m = Multigrid();
			dz0 = zeros(size(g));
			adx = upwind_kernel(cvec(obj.pmu.vx));
			ady = upwind_kernel(cvec(obj.pmu.vy));
			ad(1,:,:) = shiftdim(adx',-1);
			ad(2,:,:) = shiftdim(ady',-1);
			obj.aux.ad = ad;
			obj.aux.e  = [rvec(obj.pmu.ex); rvec(obj.pmu.ey)];
			obj.mg_m.init(obj.aux.aa,obj.aux.ad,obj.aux.e,obj.L,obj.nx,obj.nvar);
		case {'multigrid-java'}
			obj.aux.compute_aa = true;
			init_multigrid_java();
			obj.aux.compute_A = false;
		otherwise
			error('umimplemented inner-solver');
		end
		obj.aux.fstep = @obj.step_implicit;
	otherwise
		% build-in solver
		obj.init_advection_diffusion_matrix();
	end % switch obj.opt.solver
	obj.out.runtime(1) = toc();

	function init_multigrid_java()
		obj.aux.mg_j = javaObject('Multigrid_java');
		obj.aux.mg_j.set_reltol(obj.opt.inner2_tol);
		obj.aux.mg_j.set_nmaxiter(obj.opt.inner2_maxiter);
		obj.aux.dz0 = zeros(obj.nvar,prod(obj.nx));
		obj.aux.tflag = true;
		obj.jacobian_react(0,obj.aux.dz0(:),obj.aux.dz0,2,obj.aux.tflag); 
		adx = upwind_kernel(cvec(obj.pmu.vx));
		ady = upwind_kernel(cvec(obj.pmu.vy));
		ad(1,:,:) = shiftdim(adx',-1);
		ad(2,:,:) = shiftdim(ady',-1);
		ad = zeros(2,3,obj.nvar);
		for idx=1:obj.nvar
			ad(2,:,idx) = upwind_kernel(-obj.pmu.vx(idx));
			ad(1,:,idx) = upwind_kernel(-obj.pmu.vy(idx));
		end
		obj.aux.ad = -ad;
		obj.aux.e  = [rvec(obj.pmu.ey); rvec(obj.pmu.ex)];
		obj.aux.mg_j.init(obj.aux.aa,obj.aux.ad,obj.aux.e,obj.L,obj.nx);
		obj.aux.isaa = true;
	end
end

