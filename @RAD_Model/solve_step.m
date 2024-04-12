% Mon  2 May 14:18:38 CEST 2022
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
	% initial time
	% TODO allow for restart
	t = 0;
	% initial time step
	dt  = min(obj.opt.dt,obj.T);
	% output time step
	dto = min(max(dt,obj.opt.dto),obj.T);
	%tt  = (0:dt:obj.T);
	to  = (0:dto:obj.T);
	no = length(to);
	nt = no;

	obj.out.runtime     = zeros(1,no);
	obj.out.rmse        = zeros(1,no);
	obj.out.n_step      = zeros(1,no);
	obj.out.esum        = zeros(1,no);
	obj.out.dt          = zeros(1,no);
	% allocate aux output for implicit solver
%	if (strcmp(obj.opt.solver,'solve_implicit'))
		%obj.out.innersteps  = zeros(nt,1);
		%obj.out.exitflag    = zeros(nt,1);
		%obj.out.exitflag2   = {[NaN,NaN]};
		%obj.out.innersteps2 = {[NaN,NaN]};
		%obj.out.fs  = zeros(nt,1);
		%obj.out.a_line_search = [];
		%obj.out.flag2 = {};
%	end

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

	% allocate memory for output
	% value at last 3 time steps stored for error control
	zz = zeros(obj.nvar*prod(obj.nx),nc,func2str(obj.opt.compute_class));
	tt = zeros(1,nc);
	
	% for output
	zo = zeros_man(no,obj.nvar*prod(obj.nx),func2str(obj.opt.output_class));
	% initial value
	tt(1)   = t;
	zz(:,1) = obj.z0;
	zo(1,:) = obj.z0;

	% first element contains stat value, so at first iteration tdx=2
	% time step counter
	tdx = 1;
	% output counter
	odx = 1;
	timer = tic();
	esum  = 0;
	dtold = inf;
	while (t < obj.T)
		tdx = tdx+1;
		% reallocate, double storage
		if ( tdx>2 || tdx>nt)
			if (tdx > 2)
				nt    = 2*nt;
			end
	%		obj.out.stat(nt) = obj.out.stat(1);
			%tt(nt) = 0;
			%obj.out.rmse(nt) = 0;
			%obj.out.fs(nt)   = 0;
			%obj.out.a_line_search(nt) = 0;
			%obj.out.innersteps(nt) = 0;
			%obj.out.exitflag(1:nt) = 0;
			%obj.out.exitflag2{nt} = {};
			%obj.out.inntersteps2{nt} = {};
		end

		% step from t, tdx-1 to t+dt, tdx
		timeri = tic();
		%[zz(:,cc(1)),obj.out.stat(tdx)] = obj.aux.fstep(t,zz(:,cc(2)),dt,tdx,tt,zz);
		% TODO, this is inefficient, as it reorders columns, so better pass cc
%		cc_ = circshift(cc,1);
%		cc_old = cc;
%		cc = circshift(cc,+1);
%		tt_(cc_) = tt;
%		zz_(:,cc_) = zz;
		zold = zz(:,1);
		eo = obj.opt.extrapolation_order;
		eo = min(eo,tdx);
		c = (vander_1d(dt,eo)*inv(vander_1d(tt(1:(eo+1))-tt(1),eo)))';
		z = zz(:,1:(eo+1))*c;
		% TODO flag if variables are positive
		z = max(0,z);
		
		if (0)
		switch (obj.opt.extrapolation_order)
		case {0} % constant
			z = zz(:,1); 
		case {1} % linear
			%z = zz(:,1) + (dt/dtold)*(zz(:,1)-zz(:,2));
			if (tdx>2)
				c = (vander_1d(dt,1)*inv(vander_1d(tt(1:2)-tt(1),1)))';
				z = zz(:,1:2)*c;
			else
				z=zz(:,1);
			end
			% TODO flag if variables are positive
			z = max(z,0);
		case {2} % quadratic
			if (tdx>3)
				c = (vander_1d(dt,2)*inv(vander_1d(tt(1:3)-tt(1),2)))';
				z = zz(:,1:3)*c;
			else
				z=zz(:,1);
			end
			z = max(z,0);
		otherwise
			error('not yet implemented');
		end
		end

		tt = circshift(tt,1);
		zz = circshift(zz,1,2);

		attempt = 0;
		stat = [];
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
		stat.attempt=attempt;
		t            = t+dt;
		tt(1)        = t;
%		tt(cc(1))    = t;
		stat.runtime = toc(timeri);
		if (tdx>=nc)
			esum = esum+stat.rmse;
		end

		% store output
		if (t >= to(odx+1) || stat.flag(end))
			printf('Progress: %f%% Time: %f dt: %f dt0: %f \n',100*t/obj.T,t,dt,stat.dt0);
			odx = odx+1;
			% reallocate output array (this should not be required)
			if (odx>no)
				no              = 2*no;
				to(no)          = 0;
				zo(1,no)        = 0;
				obj.out.rmse(no)    = 0;
				obj.out.runtime(no) = 0;
			end
			to(odx)   = t;
			%zo(odx,:) = zz(:,cc(1));
			zo(odx,:) = zz(:,1);
			% since the first column is the initial value steps = tdx-1
			obj.out.n_step(odx) = tdx-1;
			obj.out.runtime(odx) = toc(timer);
			obj.out.rmse(odx) = stat.rmse; %rms(e);
			obj.out.stat(odx) = stat;
			obj.out.esum(odx) = esum;
			obj.out.dt(odx) = dt;
%			if (obj.opt.condest)
%				obj.out.condest(odx) = condest(obj.aux.AA);
%			end
		end
		b = obj.extract2(zz(:,1));
%		imagesc(b)
%		dt
%		t
%		stat
%pause(1)
%		stat.inner2
		if (0 ~= stat.flag(end))
			finish();
			% in case of error store last states
			obj.out.tt = tt;
			obj.out.zz = zz;
			break;
		end

		% adaptive error control, determine the time step
%dt
		dtold = dt;
		if (obj.opt.adapt_time_step && tdx>=nc)
			dt = obj.adapt_time_step(t,dt,stat.dt0,stat);
		end
%dt
%pause
	end % for tdx
	obj.out.y_final = zz(:,cc(1));

	% shrink matrices to actual size
	nt = tdx;
%	if (strcmp(obj.opt.solver,'solve_implicit'))
%	obj.out.a_line_search = obj.out.a_line_search(1:nt);
%	obj.out.exitflag      = obj.out.exitflag(1:nt);
%	obj.out.exitflag2     = obj.out.exitflag2(1:nt);
%	obj.out.fs            = obj.out.fs(1:nt);
%	obj.out.innersteps    = obj.out.innersteps(1:nt);
%	obj.out.innersteps2   = obj.out.innersteps2(1:nt);
%	obj.out.step.runtime       = obj.out.step.runtime(1:nt);
	%obj.out.t             = tt(1:nt);
%	end
	finish();

function finish()
	no                = odx
	to                = to(1:no);
	zo                = zo(1:no,:);
	obj.out.rmse      = obj.out.rmse(1:no);
	obj.out.runtime   = obj.out.runtime(1:no);
	obj.out.n_step    = obj.out.n_step(1:no);
	obj.out.dt        = obj.out.dt(1:no);
	out = obj.out;
end
end % solve_step

