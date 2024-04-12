% Thu  7 Dec 16:41:17 CET 2023
function dt = adapt_time_step(obj,t,dt_old,dt0,stat)
	%abstol  = obj.opt.outer_abstol;
	%maxiter = obj.opt.inner_maxiter;
	%rmse_div_dt2 = rms(e_div_dt2);
	% max time step satisfying the error constraint
	%dt = sqrt(abstol./rmse_div_dt2)
	dt = dt0;

	% limit at lower bound
	dt = max(dt,obj.opt.dt_min);

	% limit at upper bound
	dt = min(dt,obj.opt.dt_max);

 	% limit decrease
	dt = max(dt,dt_old/obj.opt.dt_scale);

	% limit increase
	% TODO this can be too high when the time step was reduced
	% do to a failed innersolver attemp
	dt = min(dt,dt_old*obj.opt.dt_scale);

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
			%emin = min(stat.emin(end,:));
			%minaa = stat.minaa(end);
			%a = (emin-1)/dt;
			% a*dt_min > -1
			%if (a < 0)
			%dt_max = -0.99/a;
			%dt0 = min(dt0,dt_max);
			%end
		end
	end % switch i2
	end % switch innersolver

	% limit final step to final time
	dt = min(dt,obj.T-t);

%	if (dt == 0)
%		error('adapt_time_step');
%	end
	% limit time step to hit next output time
	% dt_ = min(dt,to(todx+1)-t);
	%dt_ = dt;
end

