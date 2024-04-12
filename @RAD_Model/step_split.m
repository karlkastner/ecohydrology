% Mon  2 May 14:18:38 CEST 2022
function [z, stat] = step_split(obj,t,z,zold,dt,tt,zz)
	% react half step
	%zold = z;
	z = step_react_heun(t,0.5*dt,z,@obj.dz_dt_react);
	% advect-diffuse a full step
	z = step_advect_diffuse(t+0.5*dt,z);
	% react half step
	[z,e] = step_react_heun(t+0.5*dt,0.5*dt,z,@obj.dz_dt_react);

	% bc we step two-half steps
	rmse  = 2*rms(e);
	erel  = rmse/rms(z);
	erel0 = obj.opt.outer_reltol;
	dt0   = dt*cbrt(erel0./erel);

	stat = struct('rmse',2*rms(e),'dt0',dt0,'flag',0);

	function z = step_advect_diffuse(t,z)
		%zold = z;
		%z_C = cell(1,obj.nvar);
		%[z_C{:}] = obj.extract2(z);
		if (1 == obj.ndim)
			z_ = reshape(z,[obj.nx(1),obj.nvar]);
		else
			z_ = reshape(z,[obj.nx(1),obj.nx(2),obj.nvar]);
		end

		switch (obj.opt.inner_solver)
		case {'step_advect_diffuse_trapezoidal'}
		for idx=1:obj.nvar
			if (1 == obj.ndim)
				%z_C{idx} = ifft(obj.aux.F{idx}.*fft(z_C{idx}));
				z_(:,:,idx) = ifft(obj.aux.F{idx}.*fft(z_(:,:,idx))); %C{idx}));
			else
				%z_C{idx} = flat(ifft2(obj.aux.F{idx}.*fft2(z_C{idx})));
				z_(:,:,idx) = ifft2(obj.aux.F{idx}.*fft2(z_(:,:,idx)));
			end % if
			if (obj.opt.isreal)
				z_ = real(z_);
			end
		end
		case {'step_advect_diffuse_spectral'}
			dx = obj.dx;
			for idx=1:obj.nvar
				if (1 == obj.ndim)
				%z_C{idx} = flat(step_advect_diffuse_spectral(dt,dx, ...
				%		[obj.pmu.vx(idx),obj.pmu.vy(idx)], ...
				%		[obj.pmu.ex(idx),obj.pmu.ey(idx)], ...
				%		z_C{idx}));
				z_(:,idx) = step_advect_diffuse_spectral(dt,dx, ...
						[obj.pmu.vx(idx)], ...
						[obj.pmu.ex(idx)], ...
						z_(:,idx),obj.opt.isreal);
				else
				z_(:,:,idx) = step_advect_diffuse_spectral(dt,dx, ...
						[obj.pmu.vy(idx),obj.pmu.vx(idx)], ...
						[obj.pmu.ey(idx),obj.pmu.ex(idx)], ...
						z_(:,:,idx),obj.opt.isreal);
				end
			end
		case {'step_advect_diffuse_implicit_q_fft'}
			dx = obj.dx;
			for idx=1:obj.nvar
				if (1 == obj.ndim)
				z_(:,idx) = step_advect_diffuse_implicit_q_fft( ...
					 dt,dx,[obj.pmu.vx(idx)], ...
					    [obj.pmu.ex(idx)], ...
					    z_(:,idx),obj.opt.inner_q,obj.opt.isreal);
				%z_C{idx} = real(flat(step_advect_diffuse_implicit_q_fft( ...
				%	 dt,[obj.pmu.vy(idx),obj.pmu.vx(idx)], ...
				%	    [obj.pmu.ey(idx),obj.pmu.ex(idx)], ...
				%	    z_C{idx},obj.opt.inner_q)));
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
		%z = vertcat(z_C{:});
	end	
end % step_split

