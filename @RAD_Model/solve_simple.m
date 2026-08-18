% Thu 11 Jun 19:03:05 CEST 2026
function [to, zo, out] = solve(obj)

	dt = obj.opt.time_integration.dt_initial;
	
	timer = tic();
	
	told = 0;
	iold = 0;
	stat  = struct();
	
	nto = ceil((obj.T(2)-obj.T(1))./obj.opt.output.dt);
	
	% allocate memory
	zo = zeros(nto+1,prod(obj.nx),func2str(obj.opt.compute_class));
	to = obj.T(1) + (0:nto)*(obj.T(2)-obj.T(1))./(nto);
	% prevent overflow
	to(end)=inf;
	
	nt = ceil((obj.T(2)-obj.T(1))./dt);
	dt = (obj.T(2)-obj.T(1))./nt;
	b = obj.z0;
	t  = 0;
	odx = 1;
	dx = obj.dx;
	zo(1,:) = b(:);
	for idx=1:nt
		% splitting scheme
		% half reaction step
		b = step_react_heun(t,0.5*dt,b,@obj.dz_dt_react);
		% full diffusion step
		if (obj.ndim == 1)
			%b = (obj.aux.I - (0.5*dt*obj.p.D(1))*obj.aux.D2) \ (b + (0.5*dt*obj.p.D(1))*(obj.aux.D2*b));
			b = step_advect_diffuse_implicit_q_fft(dt,dx,[0,0],obj.p.ex{1},b,obj.opt.time_integration.q);  
		else
			b = step_advect_diffuse_implicit_q_fft(dt,dx,[0,0],[obj.p.ex{1},obj.p.ey{1}],b,obj.opt.time_integration.q);  
		end % switch ndim
		% half reaction-step
		b = step_react_heun(t,0.5*dt,b,@obj.dz_dt_react);
		b = max(b,0);
		if (t >= to(odx+1))
			printf('T = %g progress = %g%%\n',t,100*t/(obj.T(2)-obj.T(1)));
			odx = odx+1;
			to(odx) = t;
			zo(odx,:) = b(:);
		end % if t>to(odx+1)
		% advance to next time
		t = t+dt;
	end % for idx t
	odx = odx+1;
	to(odx) = t;
	zo(odx,:) = b(:);
	
	obj.out.runtime = toc(timer);
	obj.out.to = to;
	obj.out.zo = zo;
	obj.out.z_final = b;
	out = obj.out;
end % solve

