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
function [z, stat] = step_split(obj,t,z,zold,dt,tt,zz)

	% react half step
	z = step_react_heun(t,0.5*dt,z,@obj.dz_dt_react);

	sigma = obj.sigma(t,z); 

	% advect-diffuse a full step
	z = step_advect_diffuse(t+0.5*dt,z);

	% account for temporal noise
	if (~isempty(sigma))
		% TODO allow for correlation in time
		dz = sqrt(dt)*sigma.*randn(obj.nvar*prod(obj.nx),1);
		z = (z + dz);
	end

	% react half step
	[z,maxe,dt_opt] = step_react_heun(t+0.5*dt,0.5*dt,z,@obj.dz_dt_react,obj.opt.outer_abstol,obj.opt.outer_reltol);

	% conservative estimates, since the step is split in two, we assume the
	% note that the splitting scheme is of similar order of accuracy as
	% Heun's scheme, we thus use the error of the reaction part only
	% for the error estimate
	stat = struct('maxe',maxe,'dt_opt',dt_opt,'flag',0);

	function z = step_advect_diffuse(t,z)
		if (1 == obj.ndim)
			z_ = reshape(z,[obj.nx(1),obj.nvar]);
		else
			z_ = reshape(z,[obj.nx(1),obj.nx(2),obj.nvar]);
		end

		switch (obj.opt.inner_solver)
		case {'step_advect_diffuse_trapezoidal'}
			% computed in frequency space
			for idx=1:obj.nvar
				if (1 == obj.ndim)
					z_(:,:,idx) = ifft(obj.aux.F{idx}.*fft(z_(:,:,idx)));
				else
					z_(:,:,idx) = ifft2(obj.aux.F{idx}.*fft2(z_(:,:,idx)));
				end % if
				if (obj.opt.isreal)
					z_ = real(z_);
				end
			end % for idx
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
						, obj.opt.isreal);
				else
					z_(:,:,idx) = step_advect_diffuse_spectral( ...
						z_(:,:,idx), ...
						dt, ...
						[obj.pmu.vy(idx),obj.pmu.vx(idx)], ...
						[obj.pmu.ey(idx),obj.pmu.ex(idx)], ...
						obj.nx, ...
						obj.L ...
						,obj.opt.isreal);
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
					 dt,dx,[obj.pmu.vx(idx)], ...
					    [obj.pmu.ex(idx)], ...
					    z_(:,idx),obj.opt.inner_q,obj.opt.isreal);
				else
					z_(:,:,idx) = step_advect_diffuse_implicit_q_fft( ...
					 dt,dx,[obj.pmu.vx(idx),obj.pmu.vy(idx)], ...
					    [obj.pmu.ex(idx),obj.pmu.ey(idx)], ...
					    z_(:,:,idx),obj.opt.inner_q,obj.opt.isreal);
				end
			end
		otherwise
			error('unimplemented scheme');
		end
		% TODO instead of reshaping it at every time we leave it as a matrix
		z = reshape(z_,obj.nvar*prod(obj.nx),1);
	end	
end % step_split

