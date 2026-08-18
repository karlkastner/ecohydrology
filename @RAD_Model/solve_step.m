% 2021-07-06 10:56:59.150856622 +0200
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
function solve_step(obj)
	% initial time
	t = obj.T(1);
 	% initial state
	z = obj.z0;
	% initial time step
	dt  = min(obj.opt.time_integration.dt_initial,obj.T(2));

	% for output at relative change
	% note that zo_last has to be stored in the same precision as z and not
	% in the output precision, esp when the type of zo is this can otherwise result in overflow
	obj.aux.zo_last = obj.z0;
	obj.aux.rms_zo_last = rms(obj.aux.zo_last);
	obj.aux.rms_delta_zo = 0;

	% first element contains the start value, so the values at t=dt
	% is written to tdx=2
	% time step counter
	obj.aux.tdx = 1;
	% output counter
	obj.aux.odx = 1;
	% event counter
	obj.aux.edx    = 1;
	% timer
	obj.aux.timer = tic();
	while (t < obj.T(2))
		% step preprocess
		[z] = obj.step_preprocess(t,dt,z);

		timeri = tic();
		obj.aux.tdx  = obj.aux.tdx+1;
		obj.aux.zold = z;
		obj.aux.told = t; 
		obj.aux.t    = t;

		if (obj.advance_to_next_event(t,z))
			dt = max(dt,obj.aux.tevent(obj.aux.edx+1)-t);
		end
	
		% limit final step to final time
		dt = min(dt,obj.T(2)-t);

		if (t+dt >= obj.aux.tevent(obj.aux.edx+1))
			% step only until next event
			% nb: the min is not necessary, as always smaller dt
			dt_ = min(dt,obj.aux.tevent(obj.aux.edx+1)-t);
			% increase event counter after current time step
			dedx = 1;
		else
			dedx = 0;
			dt_ = dt;
		end

		obj.aux.step.t = t;
		obj.aux.step.n_attempt = 0;
		obj.aux.step.n_error_tolerance_exceeded = 0;
		obj.aux.step.n_solver_failed = 0;
		obj.aux.step.n_neg = 0;
		obj.aux.step.n_iter = 0;
		obj.aux.step.n_liter = 0;
		% step, if implicit solver fails, reduce the time step
		while (1)
			obj.aux.step.n_attempt = obj.aux.step.n_attempt+1;
			% store dt as aux variable, this is needed for contant rate	
			obj.aux.dt = dt_;
			[z,obj.aux.stat] = obj.aux.fstep(t,z,dt_);
			obj.aux.stat.dt = dt_;
			if (isfield(obj.aux.stat,'iter'))
				%citer = citer + stat.iter;
				obj.aux.step.n_iter = obj.aux.step.n_iter + obj.aux.stat.iter;
			end
			if (isfield(obj.aux.stat,'linear_solver'))
				obj.aux.step.n_liter = obj.aux.step.n_liter + sum([obj.aux.stat.linear_solver.iter]);
			end
			% TODO no magic numbers and ranges
			%if (min(z(1:2*end/3)) < -1e-4)
			if (min(z) < -1e-4)
				printf('Progress %f%% Time %f Step %d Negative value, reducing time step to %g\n' ...
					, 100*t/obj.T(2), t,obj.aux.tdx,dt_);
				z=obj.aux.zold;
				dt_ = dt_/sqrt(2);
				dedx = 0;
				obj.aux.step.n_neg = obj.aux.step.n_neg +1;
			elseif (obj.aux.stat.flag(end))
				obj.aux.step.n_solver_failed = obj.aux.step.n_solver_failed + 1;
				% solver failed, reduce time step and try again
				% note : this should almost never happen
				dt_ = dt_/sqrt(2);
				dt  = dt_;
				dedx=0;
				% reset starting point
				z = obj.aux.zold;
				printf('Progress %f%% Time %f Step %d Solver failed, reducing time step to %g\n' ...
					, 100*t/obj.T(2),t,obj.aux.tdx,dt_);
			elseif (obj.opt.time_integration.adapt_step_length ...
				&& obj.aux.stat.dt_opt < 0.8*dt_ ...	
				&& dt_ > obj.opt.time_integration.dt_min ...
				)
				obj.aux.step.n_error_tolerance_exceeded = obj.aux.step.n_error_tolerance_exceeded + 1;
				% if the error is too large, decrease the time step and retry
				dt_ = max(obj.opt.time_integration.dt_min,obj.aux.stat.dt_opt);
				dedx=0;
				% interpolate as initial guess
				z = obj.aux.zold; % + obj.aux.stat.dt_opt/dt_*(z-obj.aux.zold);
				printf('Progress %f%% Time %f Step %d Error tolerance exceeded, reducing time step to %g\n' ...
					, 100*t/obj.T(2),t,obj.aux.tdx,dt_);
			else
				break;
			end
			if (obj.aux.step.n_attempt > obj.opt.time_integration.max_attempts || dt_==0)
				disp('Too many attempts, terminating');
				obj.aux.stat.runtime = toc(timeri);
				obj.step_postprocess(z);
				obj.store(t,z);
				obj.finish(z);
				obj.out.z_final = NaN(obj.nvar*prod(obj.nx),1);
				obj.out.zold = obj.aux.zold;
				return;
			end
		end % while 1
		obj.aux.step.runtime = toc(timeri);

		% store step statistics
		obj.step_postprocess(z);

		% step in time
		t            = t+dt_;

		% difference of current state to last stored state
		if (isfinite(obj.opt.rms_delta_zo_rel_max))
			obj.aux.rms_delta_zo = rms(z - obj.aux.zo_last);
		end

		% store output of time since last storing exceeds dto
		% or change since last time step exceeds delta_zo_rel_max
		if (	(obj.aux.odx+1 < length(obj.out.to) && (t >= obj.out.to(obj.aux.odx)+obj.opt.output.dt)) ...
		    || (obj.aux.rms_delta_zo >= obj.opt.rms_delta_zo_rel_max*obj.aux.rms_zo_last) ...
		    ... % || (obj.output_event(t,z,dt)) ...
		   )
		    % ... || obj.aux.stat.flag(end)  ...
			printf('Progress %f%% Time %f Step %d ||dzo||/||zo||: %g dt: %f dt_opt: %f dt_max_prec: %f\n' ...
					, 100*t/obj.T(2) ...
					, t ...
					, obj.aux.tdx ...
					, obj.aux.rms_delta_zo/obj.aux.rms_zo_last ...
					, dt ...
					, obj.aux.stat.dt_opt ...
					, obj.aux.dt_max_prec ...
			      );
			obj.store(t,z);
		end

		if (0 ~= obj.aux.stat.flag(end))
			% note this should neve happen, since errors are captured above
			obj.finish(z);
			break;
		end

		% adaptive error control, determine the time step
		if (obj.opt.time_integration.adapt_step_length)
			dt = obj.adapt_time_step(t,dt);
		else
			% reset to the initial time step
			% due to retries the time step might have been adapted
			% temporarily
			dt = obj.opt.time_integration.dt_initial; 	
			% max(min(dt,obj.opt.time_integration.dt_max),obj.opt.time_integration.dt_min);
		end
		if (obj.opt.time_integration.adapt_step_for_stability)
			dt = min(dt,obj.aux.stability_factor*obj.aux.dt_max);
		end

		if (dedx)
			% limit time step at events
			% (jump of source terms in time)
			dt = min(dt,obj.opt.time_integration.dt_event);
			obj.aux.edx = obj.aux.edx+dedx;
			dedx = 0;
		end
	end % while t < T(end)

	% write final step
	if (obj.out.to(obj.aux.odx)<obj.T(2))
		obj.store(t,z);
	end

	obj.finish(z);

end % solve_step

