% Wed 13 Mar 15:32:05 CET 2024
	param = struct();
	param.L  = 128*[1,1];
	param.nx = 128*[1,1];
	param.T  = 150;

solver_C = {'euler_forward' ...
	    ,'solve_split' ...
	    ,'solve_split' ...
	    ,'step_integrating_factor' ...
	    ,'step_react_advect_diffuse_erk' ...
... %	    ,'solve_implicit' ...
	   };

inner_solver_C = {        '' ...
			, 'step_advect_diffuse_spectral' ...
			, 'step_advect_diffuse_implicit_q_fft' ...
			, '' ...
			, '' ... 
... %			, 'gauss-newton' ...
		      };


np = length(solver_C)
for idx=1:np

	param.opt.dto    = 1; 
        param.opt.dt     = 0.25/4;
	param.opt.adapt_time_step = 0; %adapt_time_step(idx);
	param.opt.rng = 1;
	param.opt.inner_solver = inner_solver_C{idx};
	param.opt.solver = solver_C{idx}; %'solve_split';	
	param.opt.path_str = './mat/';
	param.opt.output_class = @single;
	param.opt.compute_class = @double;

	rad = Ginzburg_Landau(param);
	
	[t,z] = rad.run();

	subplot(2,np,idx);
	cla();
	u = rad.extract2(z(1,:));
	%imagesc(real(reshape(z(end,:),rad.nx)))
	imagesc(real(double(u)));
	colorbar
	title(solver_C{idx},inner_solver_C{idx})

	subplot(2,np,np+idx);
	cla();
	u = rad.extract2(z(end,:));
	%imagesc(real(reshape(z(end,:),rad.nx)))
	imagesc(real(double(u)));
	colorbar

	

end

