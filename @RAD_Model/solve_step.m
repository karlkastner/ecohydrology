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
function [to,zo,out] = solve_step(obj)
	T = obj.T;
	if (isscalar(T))
		T = [0,T];
	end

	% initial time
	t = T(1);
	% initial time step
	dt  = min(obj.opt.dt,T(2));
	% output time step
	% note that the output time step will vary when output is set to be
	% written when the state variable exceeds the maximum relative change
	if (~isempty(obj.opt.dto))
		dto = min(max(dt,obj.opt.dto),T(2)-T(1));
	else
		dto = T(2)-T(1);
	end
	to  = (T(1):dto:T(2));
	no  = length(to);
	nt  = no;

	obj.out.runtime     = zeros(1,no);
	obj.out.rmse        = zeros(1,no);
	obj.out.n_step      = zeros(1,no);
	obj.out.esum        = zeros(1,no);
	obj.out.dt          = zeros(1,no);

	% circling counter for adaptive error control
	if (strcmp(obj.opt.solver,'solve_implicit'))
	switch (obj.opt.inner_q)
	case {1}
		nc  = 3;
	case {0.5}
		nc  = 4;
	otherwise
		nc = 1;
	end
	else
		nc = 1;
	end
	cc  = (1:nc);

	% allocate memory for current state variable, overwritten every nc time steps
	% keep last nc time steps stored for error control
	zz = zeros(obj.nvar*prod(obj.nx),nc,func2str(obj.opt.compute_class));
	tt = zeros(1,nc);
	
	% state variable at output times
	zo = zeros_man(no,obj.nvar*prod(obj.nx),func2str(obj.opt.output_class));

	% initial value
	tt(1)   = t;
	zz(:,1) = obj.z0;
	zo(1,:) = obj.z0;
	% note that zo_last has to be stored in the same precision as z and not zo
	% half can lead to overflow
	zo_last = obj.z0;
	%(odx,:)');
	rms_zo_last = rms(zo(1,:));

	% first element contains stat value, so at first iteration tdx=2
	% time step counter
	tdx = 1;
	% output counter
	odx = 1;
	timer = tic();
	esum  = 0;
	dtold = inf;
	while (t < T(2))
		timeri = tic();
		tdx  = tdx+1;
		zold = zz(:,1);

		% step from t, tdx-1 to t+dt, tdx
		%[zz(:,cc(1)),obj.out.stat(tdx)] = obj.aux.fstep(t,zz(:,cc(2)),dt,tdx,tt,zz);
		% TODO, this is inefficient, as it reorders columns, so better pass cc
%		cc_ = circshift(cc,1);
%		cc_old = cc;
%		cc = circshift(cc,+1);
%		tt_(cc_) = tt;
%		zz_(:,cc_) = zz;
		
		% initial guess for implicit solverm, polynomial extrapolation
		eo = obj.opt.extrapolation_order;
		eo = min(eo,tdx);
		c = (vander_1d(dt,eo)*inv(vander_1d(tt(1:(eo+1))-tt(1),eo)))';
		z = zz(:,1:(eo+1))*c;
		% TODO flag if variables are positive
		z = max(0,z);
		
		tt = circshift(tt,1);
		zz = circshift(zz,1,2);

		attempt = 0;
		stat = [];
		% step, if implicit solver fails, reduce the time step
		while (1)
			attempt = attempt+1;
			[zz(:,1),stat_] = obj.aux.fstep(t,z,zold,dt,tt,zz);
			stat_.dt = dt;
			stat_.attempt_list = stat;
			stat = stat_;
			if (~stat.flag(end))
				break;
			else
				dt = dt/sqrt(2);
				disp(sprintf('Solver failed, reducing time step to %g\n',dt));
			end
			z = zz(:,1);
		end
		stat.attempt = attempt;
		t            = t+dt;
		tt(1)        = t;
		stat.runtime = toc(timeri);
		if (tdx>=nc)
			esum = esum+stat.maxe;
		end

		% difference of current state to last stored state
		rms_delta_zo = rms(zz(:,1) - zo_last);

		% store output of time since last storing exceeds dto
		% or change since last time step exceeds delta_zo_rel_max
		if (	(odx+1 < length(to) && t >= to(odx)+obj.opt.dto) ...
			|| stat.flag(end)  ...
			|| (rms_delta_zo >= obj.opt.rms_delta_zo_rel_max*rms_zo_last) ...
		   )
			printf('Progress: %f%% Time: %f ||dzo||/||zo||: %g dt: %f dt_opt: %f \n', ...
					100*t/obj.T, ...
					t, ...
					rms_delta_zo/rms_zo_last, ...
					dt, ...
					stat.dt_opt ...
			      );
			store();
		end

		if (0 ~= stat.flag(end))
			finish();
			% in case of error store last states
			obj.out.tt = tt;
			obj.out.zz = zz;
			break;
		end

		% adaptive error control, determine the time step
		dtold = dt;
		if (obj.opt.adapt_time_step && tdx>=nc)
			dt = obj.adapt_time_step(t,dt,stat.dt_opt,stat);
		end
	end % for tdx
	% write last step
	if (to(odx)<T(2))
		store();
	end

	obj.out.y_final = zz(:,cc(1));

	finish();

	function store()
			% move to next output slot
			odx = odx+1;
			% reallocate output array
			% only takes effect when output is variable
			% with respect to relative change
			if (odx>no)
				% double vector length
				no                  = 2*no;
				% extend vectors by filling zeros
				% to has to be filled with NaN or inf to avoid storing every time step
				to(end+1:no)        = NaN;
				zo(1,no)            = 0;
				obj.out.runtime(no) = 0;
				obj.out.n_step(no)  = 0;
				obj.out.runtime(no) = 0;
				obj.out.rmse(no)    = 0;
				oub.out.dt(no)      = 0;
			end
			to(odx)     = t;
			zo_last     = zz(:,1);
			rms_zo_last = rms(zo_last);
			zo(odx,:)   = zz(:,1);
			% since the first column is the initial value steps = tdx-1
			obj.out.n_step(odx)  = tdx-1;
			obj.out.runtime(odx) = toc(timer);
			obj.out.stat(odx)    = stat;
			obj.out.esum(odx)    = esum;
			obj.out.dt(odx)      = dt;
	end % store

	function finish()
		% shrink matrices to size containing values
		no                = odx;
		to                = to(1:no);
		zo                = zo(1:no,:);
		obj.out.rmse      = obj.out.rmse(1:no);
		obj.out.runtime   = obj.out.runtime(1:no);
		obj.out.n_step    = obj.out.n_step(1:no);
		obj.out.dt        = obj.out.dt(1:no);
		out               = obj.out;
	end
end % solve_step

