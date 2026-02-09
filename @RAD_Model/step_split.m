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
function [z, stat] = step_split(obj,t,z,dt)
	%maxe_ = 0;
	%idmax_ = 0;
	%dt_opt_ = inf;
	%limiting_part = 0;

	dz_dt0 = obj.dz_dt(t,z); 

	% react half step
	z = step_react_heun(t,0.5*dt,z,@obj.dz_dt_react);

	if (obj.opt.temporal_noise)
		sigma = obj.sigma(t,z); 
	end

	% advect-diffuse a full step
	z = step_advect_diffuse(t+0.5*dt,z);

	% account for temporal noise
	%if (~isempty(sigma))
	if (obj.opt.temporal_noise)
		% TODO allow for correlation in time
		dz = sqrt(dt)*sigma.*randn(obj.nvar*prod(obj.nx),1);
		z = (z + dz);
	end

	% react half step
	z = step_react_heun(t+0.5*dt,0.5*dt,z,@obj.dz_dt_react);
	%[z,maxe,dt_opt,idmax] = step_react_heun(t+0.5*dt,0.5*dt,z,@obj.dz_dt_react,obj.opt.outer_abstol,obj.opt.outer_reltol);

	if (nargout()>1)
	d = (obj.dz_dt(t+dt,z) - dz_dt0);
	[maxd,idmax] = max(d);
	maxe = 0.5*dt*maxd;
	% tol = obj.opt.outer_abstol + obj.opt.outer_reltol*rms(z,'all');
	tol = obj.time_integration_tolerance(z);
	dt_opt = dt*sqrt(tol/maxe);
	% error estimate
	%[maxe_,idmax_] = max(abs(obj.aux.zero_inertia.dh_dt(t+dt,z_(:,:,idx))-dh_dt0));
	%maxe_ = 0.5*dt*maxe_;
	%dt_opt_ = dt*sqrt(tol/maxe_);

	% conservative estimates, since the step is split in two, we assume the
	% note that the splitting scheme is of similar order of accuracy as
	% Heun's scheme, we thus use the error of the reaction part only
	% for the error estimate
	%if (maxe_>maxe)
	%	idmax = idmax_+(obj.nvar-1)*prod(obj.nx);
	%	maxe = maxe_;
	%	limiting_part = 1;
	%end
	stat = struct('maxe',maxe,'dt_opt',dt_opt,'flag',0,'idmax',idmax);
	%,'limiting_part',limiting_part);
	end

	function z = step_advect_diffuse(t,z)
		if (1 == obj.ndim)
			z_ = reshape(z,[obj.nx(1),obj.nvar]);
		else
			z_ = reshape(z,[obj.nx(1),obj.nx(2),obj.nvar]);
		end

		switch (obj.opt.time_integration.scheme)
		case {'step_advect_diffuse_aid'}
			for idx=1:obj.nvar
				if (idx < obj.nvar | ~isfield(obj.opt,'nonlinear_flow') | ~obj.opt.nonlinear_flow)
				if (1 == obj.ndim)
					% for 1D aid is just regular step
					z_(:,idx) = step_advect_diffuse_implicit_q_fft( ...
					 dt,dx,[obj.pmu.vx(idx)], ...
					    [obj.pmu.ex(idx)], ...
					    z_(:,idx),obj.opt.time_integration.q,obj.opt.isreal);
				else
					z_(:,:,idx) = step_advect_diffuse_aid_cc( ...
					        z_(:,:,idx), ...
						dt, ...
						obj.aux.Ix, ...
						obj.aux.Ax1{idx}, ...
						obj.aux.Iy, ...
						obj.aux.Ay1{idx} ...
						);
				end
				else
					if (1 == obj.ndim)
						error('TODO implement (makes no sense to use aid)');
					else
						%dh_dt0 = obj.aux.zero_inertia.dh_dt(t,z_(:,:,idx));
						[z_(:,:,idx),~,~,~,~,stat.iter] = obj.aux.zero_inertia.step(t,z_(:,:,idx),dt); 
					end
				end
			end
		case {'step_advect_diffuse_spectral'}
			% analytic advection diffusion
			dx = obj.dx;
			for idx=1:obj.nvar
				if (1 == obj.ndim)
					z_(:,idx) = step_advect_diffuse_spectral( ...
						  z_(:,idx) ...
						, dt ...
						, [obj.pmu.vx(idx)] ...
						, [obj.pmu.ex(idx)] ...
						, obj.nx ...
						, obj.L ...
						);
				else
					z_(:,:,idx) = step_advect_diffuse_spectral( ...
						z_(:,:,idx), ...
						dt, ...
						[obj.pmu.vy(idx),obj.pmu.vx(idx)], ...
						[obj.pmu.ey(idx),obj.pmu.ex(idx)], ...
						obj.nx, ...
						obj.L ...
						);
				end
			end % for idx
		case {'step_advect_diffuse_implicit_q_fft'}
			% q = 0 explicit euler
			% q = 0.5 trapezoidal (implicit)
			% q = 1  implicit euler
			dx = obj.dx;
			for idx=1:obj.nvar
				if (1 == obj.ndim)
					z_(:,idx) = step_advect_diffuse_implicit_q_fft( ...
					 dt, ...
					 dx, ...
					   [obj.pmu.vx(idx)], ...
					    [obj.pmu.ex(idx)], ...
					    z_(:,idx), ...
					    obj.opt.time_integration.q ...
					   );
				else
					z_(:,:,idx) = step_advect_diffuse_implicit_q_fft( ...
					 	dt, ...
						dx, ...
						[obj.pmu.vx(idx),obj.pmu.vy(idx)], ...
					    [obj.pmu.ex(idx),obj.pmu.ey(idx)], ...
					    z_(:,:,idx), ...
					    obj.opt.time_integration.q ...
					);
				end
			end
		otherwise
			error('unimplemented scheme');
		end
		% TODO instead of reshaping it at every time we leave it as a matrix
		z = reshape(z_,obj.nvar*prod(obj.nx),1);
	end	
end % step_split

