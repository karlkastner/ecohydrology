% Tue 12 Mar 16:32:01 CET 2024
%
%	u = L u + N(u)
%	v = exp(-L*t)*u
%	u = exp(+L*t)*v
%	dv/dt = -L*exp(-L*t)*u + exp(-L*t)*du/dt
%	      = -L*exp(-L*t)*u + exp(-L*t)*(L u + N(u))
%	      = exp(-L*t)*N(u)
%             = exp(-L*t)*N(exp(+L*t)*v)
% explicit:
%	v_1 = v_0 + dt*exp(-L*t)*N(exp(+L*t)*v_0)
% implicit:
%	v_1 = 
%	u_1 = u_0 
function [u,stat] = step_integrating_factor(obj,t,u,dt,tt,zz)
	isreal_ = obj.opt.isreal;

	[fx, fy] = fourier_axis_2d(obj.L,obj.nx);
	fmax =  0.5*obj.nx./obj.L;
	p = 5/6;
	p = 1/2;
	ax = abs(fx)>=fmax(1)*(p);
	ay = abs(fy)>=fmax(2)*(p);

	%a = zeros(2,3);
	q = 0.5;
	%linfun = @(dt,dx,a,e,z0) step_advect_diffuse_implicit_q_fft(dt,dx,a,e,z0,q);
	linfun = @(dt,dx,a,e,z) step_advect_diffuse_spectral(dt,dx,a,e,z,obj.opt.isreal);

	z_C = cell(1,obj.nvar);
	% intial value 
	% v = exp(-L*tau)*u = u, tau = 0
	v = step_react_heun(0,dt,u,@vdot);
	%v = step_rk4(0,dt,u,@vdot);
	% determine u
	% u = exp(L*tau)*v, tau = dt
	[z_C{:}] = obj.extract2(v);
	for idx=1:obj.nvar
		z_C{idx} = (linfun(dt,obj.dx,[obj.pmu.vx(idx),obj.pmu.vy(idx)], ...
					     [obj.pmu.ex(idx),obj.pmu.ey(idx)],z_C{idx}));
		% de-alias
		f = fft2(z_C{idx});
		f(ax,:) = 0;
		f(:,ay) = 0;
		z_C{idx} = flat(ifft2(f));
	end
	u = vertcat(z_C{:});
	if (isreal_)
		u = real(u);
	end

	%u = max(u,0);
	% TODO error estimate
	e    = NaN;
	dt0  = dt;
	stat = struct('rmse',2*rms(e),'dt0',dt0,'flag',0);

% time derivative with integrating factor
function dv_dt = vdot(t,v)
	% solve forward in time
	[z_C{:}] = obj.extract2(v);
	for idx=1:obj.nvar
		% u  = exp(L*t)*v;
		z_C{idx} = flat(linfun(+t,obj.dx,[obj.pmu.vx(idx),obj.pmu.vy(idx)],[obj.pmu.ex(idx),obj.pmu.ey(idx)],z_C{idx}));
	end
	u = vertcat(z_C{:});
	% nonlinear (reaction) part
	u = obj.dz_dt_react(t,u);
	% solve backward in time
	zz = {};
	[z_C{:}] = obj.extract2(u);
	for idx=1:obj.nvar
		% v = exp(-L*dt)*u
		z_C{idx} = flat(linfun(-t,obj.dx,[obj.pmu.vx(idx),obj.pmu.vy(idx)],[obj.pmu.ex(idx),obj.pmu.ey(idx)],z_C{idx}));
	end
	dv_dt = vertcat(z_C{:});
end

end

