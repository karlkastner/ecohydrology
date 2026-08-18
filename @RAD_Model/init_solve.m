% 2024-03-24 11:22:25.248722605 +0100
% Karl Kästner, Berlin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
% TODO remove deprecated multigrid code
function init_solve(obj)
	if (isfield(obj.opt,'tevent'))
		obj.aux.tevent = obj.opt.tevent;
	else
		obj.aux.tevent = [0,inf];
	end
	obj.aux.dt_max_prec = NaN;
	obj.aux.max_eig_J = sqrt(eps);

	% output times
	% note that the output time step will vary when output is set to be
	% written when the state variable exceeds the maximum relative change
	if (~isempty(obj.opt.output.dt))
		%dto = min(max(dt,obj.opt.dto),obj.T(2)-obj.T(1));
		dto = min(obj.opt.output.dt,obj.T(2)-obj.T(1));
	else
		dto = obj.T(2)-obj.T(1);
	end
	% allocate output
	obj.out = struct();
	obj.out.to  = (obj.T(1):dto:obj.T(2))';
	if (isfield(obj.opt.output,'no'))
	no  = obj.opt.output.no;
	else
	no  = length(obj.out.to);
	end
	%nt  = no;

	% allocate memory for state variable at output times
	obj.out.zo = zeros(no,obj.nvar*prod(obj.nx),func2str(obj.opt.output.class));

	% first value in output is the initial value
	obj.out.zo(1,:) = obj.z0;

	% these counters are per step and have small values
	obj.out.n_attempt = zeros(no,1,'uint16');
	obj.out.n_error_tolerance_exceeded = zeros(no,1,'uint16');
	obj.out.n_solver_failed = zeros(no,1,'uint16');
	obj.out.n_neg     = zeros(no,1,'uint16');
	obj.out.n_iter    = zeros(no,1,'uint16');
	obj.out.n_liter   = zeros(no,1,'uint16');
	obj.out.n_step    = zeros(no,1,'uint16');
	obj.out.runtime   = zeros(no,1);
	obj.out.esum      = zeros(no,1,func2str(obj.opt.compute_class));

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

	% TODO this should be moved to the Rietkerk class
	if (isfield(obj.opt,'nonlinear_flow') && obj.opt.nonlinear_flow)
		switch (obj.ndim)
		case {1}
		obj.aux.zero_inertia = SWE_Zero_Inertia_1d(obj.opt.zero_inertia);
		case {2}
		obj.aux.zero_inertia = SWE_Zero_Inertia_2d(obj.opt.zero_inertia);
		end % switch ndim
		zb = obj.p.zb;
		if (isscalar(zb))
			nx = obj.nx;
			if (isscalar(nx))
				nx =[nx,1];
			end
			zb = repmat(zb,nx+2);
		end % isscalar zb
		obj.aux.zero_inertia.eps = 1e-6*obj.opt.zero_inertia.input_factor;
        	obj.aux.zero_inertia.zb   = zb;
		obj.aux.zero_inertia.returnmat = true;
   		%obj.aux.zero_inertia.C    = obj.pmu.Chezy;
   		obj.aux.zero_inertia.lcd  = obj.pmu.lcd;
       		obj.aux.zero_inertia.L    = obj.L;
	        obj.aux.zero_inertia.n    = obj.nx;
	        obj.aux.zero_inertia.boundary_condition  = obj.boundary_condition;

       	end % if nonlinear flow

	obj.aux.q = 0;
	switch (solver_str)
	case {'solve_step'}
	switch (obj.opt.time_integration.scheme)
	case {'aid','step_aid'}
		obj.aux.q = 0.5;
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_aid;
	case {'step_integrating_factor'}
		obj.aux.q = 0;
		obj.aux.fstep = @obj.step_integrating_factor;
	case {'step_react_advect_diffuse_erk'}
		obj.aux.q=0;
		obj.aux.fstep = @obj.step_react_advect_diffuse_erk;
	case {'step_euler_forward'}
		obj.aux.q = 0;
		% TODO mobe into butcher table
		obj.aux.stability_factor = 1*obj.opt.time_integration.stability_safety_factor;
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @step_euler_forward;
		obj.aux.butcher_table = butcher_table('forward');
	case {'implicit-euler','euler-implicit','step_euler_implicit'}
		obj.aux.q = 1;
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_euler_implicit;
		obj.aux.butcher_table = butcher_table('backward');
	case {'euler-implicit-2x'}
		obj.aux.q = 1;
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_euler_implicit_2x;
	case {'trapezoidal-2x'}
		obj.aux.q = 0.5;
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_trapezoidal_2x;
	case {'trapezoidal','step_trapezoidal'}
		obj.aux.q = 0.5;
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_trapezoidal;
		obj.aux.butcher_table = butcher_table('trapezoidal');
	case {'trapezoidal_explicit','step_trapezoidal_explicit'}
		obj.aux.q = 0.5;
		% TODO mobe into butcher table
		obj.aux.stability_factor = 2*obj.opt.time_integration.stability_safety_factor;
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_trapezoidal;
		obj.aux.butcher_table = butcher_table('trapezoidal');
	case {'trbdf2'}
		obj.aux.q = (2-sqrt(2))/2; 
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_trbdf2;
		obj.aux.butcher_table = butcher_table('trbdf2');
	case {'midpoint','step_midpoint'}
		obj.aux.q = 0.5;
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_midpoint;
		obj.aux.butcher_table = butcher_table('midpoint');
	case {'midpoint_explicit','step_midpoint_explicit'}
		obj.aux.q = 0.5;
		% TODO mobe into butcher table
		obj.aux.stability_factor = 2*obj.opt.time_integration.stability_safety_factor;
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_midpoint_explicit;
		obj.aux.butcher_table = butcher_table('midpoint');
	case {'sdirk23'}
		obj.aux.q=(3+sqrt(3))/6;
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_sdirk23;
		obj.aux.butcher_table = butcher_table('sdirk23');
	case {'sdirk2'}
		obj.aux.q = (2-sqrt(2))/2;
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_sdirk2;
		obj.aux.butcher_table = butcher_table('sdirk2');
	case {'euler-double','step_euler_double'}
		obj.aux.q = 1;
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_euler_double;
	case {'step_advect_diffuse_q'}
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_split;
	case {'step_advect_diffuse_aid'}
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_split;
	case {'step_advect_diffuse_spectral'}
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_split;
		% nothing to do
	case {'step_advect_diffuse_implicit_q_fft'}
		obj.init_advection_diffusion_matrix();
		obj.aux.fstep = @obj.step_split;
		% nothing to do
	%	otherwise
	%	error(sprintf('unimplemented inner solver %s',obj.opt.inner_solver));
	%	end
	otherwise
		error('Unimplemented time_stepper')
	end
	otherwise
		% build-in solver
		obj.init_advection_diffusion_matrix();
	end % switch obj.opt.solver
	obj.out.runtime(1) = toc();
end

