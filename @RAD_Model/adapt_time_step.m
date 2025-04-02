% Thu  7 Dec 16:41:17 CET 2023
function dt = adapt_time_step(obj,t,dt_old,dt0,stat)
	% max time step satisfying the error constraint
	dt = dt0;

 	% limit decrease
	dt = max(dt,dt_old*obj.opt.dt_min_scale_down);

	% limit increase
	dt = min(dt,dt_old*obj.opt.dt_max_scale_up);

	% limit at lower bound
	dt = max(dt,obj.opt.dt_min);

	% limit at upper bound
	dt = min(dt,obj.opt.dt_max);

	% avoid too many inner steps
	switch (obj.opt.inner_solver)
	case {'newton'}
	switch (obj.opt.innersolver2)
	case {'multigrid-java'}
		iter_max = obj.opt.inner2_maxiter_;
		if (~isempty(iter_max))
			% the last one is the hardest always
			iter = stat.inner2.iter(end);
			dt = min(dt,dt_old*iter_max/iter);
		end
	end % switch i2
	end % switch innersolver

	% limit final step to final time
	dt = min(dt,obj.T-t);
end

