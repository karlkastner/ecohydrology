% Mon  2 May 14:18:38 CEST 2022
%
%% evolve the Rietker-PDE in time using a splitting scheme
%
function [t,yy,fallback] = solve_split(obj,to,y0)
	dt = obj.dt;
	nto = length(to);
	nt  = round((to(end)-to(1))/dt)+1;
	t   = linspace(to(1),to(end),nt)';
	fallback = false(nto,1);
	yy  = zeros(nto,numel(y0));
	prefactor = true;
	resn = [];
	bbar = [];
	% prepare time stepping, diffusin part
	switch (obj.opt.diffusion_scheme)
	case('aid')
		% constant coefficients
		[Dx,Dy,Dxx,Dyy] = setup_matrix();
		Ax = Dx+Dxx;
		Ay = Dy+Dyy;

		if (prefactor)
		% TODO, in the test, chol was not faster or slower than LU, but required less memory
		% note : solving the system with lu-factorization takes appr.
		% as much times as 20 iterations: y = a*y+b, with a and b vectors, without setup time for a and b
		% (I-Ax) and (I-Bx) commute, so no need to alternate at odd even steps
		% Lc = chol(I-dt*Ax);
		% Lct = Lc';
		[lx,ux] = lu(I-dt*Ax);
		if (obj.ndim>1)
		[ly,uy] = lu(I-dt*Ay);
		end
		end
	case ('fdm')
		% TODO, this is only 1D
		c  = obj.dz_dt_coefficient(0,cvec(y0));
		nn = prod(obj.n);
		c3 = c{:,3};
		I  = speye(3*nn);
		D2 = [c3{1}(1)*obj.D2, obj.Z, obj.Z;
		     obj.Z, c3{2}(1)*obj.D2, obj.Z;
		     obj.Z, obj.Z, c3{3}(1)*obj.D2];
		if (obj.ndim>1)
			error('not yet implmented');
		end
	case ('spectral')
		% nothing to do, no matrices required
	case ('fdm-positive')
		% unit impulse
		z0 = zeros(prod(obj.n),1);
		z0(1,1) = 1;
		% set up the diffusion matrices
		% convolution matrices (impulse response)
		dx = obj.dx;
		gb = diffuse_trapezoidal(dt,z0,obj.n,dx,obj.p.eb);
		gw = diffuse_trapezoidal(dt,z0,obj.n,dx,obj.p.ew);
		k  = 10;

		gh = diffuse_trapezoidal(dt/k,z0,obj.n,dx,obj.p.eh);
		%gh = diffuse_euler_implicit(dt/k,z0,obj.n,dx,obj.p.eh);
	
		if (2 == obj.ndim)
			gb = reshape(gb,obj.n);
			gw = reshape(gw,obj.n);
			gh = reshape(gh,obj.n);
		end
		% ft of convolution matrices
		obj.aux.fgb = fft2(gb);
		obj.aux.fgw = fft2(gw);
		obj.aux.fgh = fft2(gh).^k;
	case {'advection-diffusion'}
		% unit impulse
		z0 = zeros(prod(obj.n),1);
		z0(1,1) = 1;
		% set up the diffusion matrices
		% convolution matrices (impulse response)
		dx = obj.dx;
		dt_max = 0.5*sum(dx.^2./obj.p.eb);
		kb = ceil(dt/dt_max);
		gb = step_diffuse_trapezoidal(dt/kb,z0,obj.n,dx,obj.p.eb);
		dt_max = 0.5*sum(dx.^2./obj.p.ew);
		kw = ceil(dt/dt_max);
		gw = step_diffuse_trapezoidal(dt/kw,z0,obj.n,dx,obj.p.ew);
		dt_max = 0.5*sum(dx.^2./obj.p.eh);
		k1  = ceil(dt/dt_max);
		k2  = ceil(max(abs(obj.p.vh)*dt./dx));
		kh   = max(k1,k2);
		kh   = max(1,kh);
%		k  = max(2,k)
		%gh = diffuse_trapezoidal(dt/k,z0,obj.n,dx,obj.p.eh);
		gh = step_advection_diffusion_trapezoidal(dt/kh,obj.dx,obj.n,z0,obj.p.vh,obj.p.eh);
		%gh = step_advection_diffusion_euler_implicit(dt/k,obj.dx,obj.n,z0,obj.p.vh,obj.p.eh);
		%gh = diffuse_euler_implicit(dt/k,z0,obj.n,dx,obj.p.eh);
	
		if (2 == obj.ndim)
			gb = reshape(gb,obj.n);
			gw = reshape(gw,obj.n);
			gh = reshape(gh,obj.n);
		end
		% ft of convolution matrices
		obj.aux.fgb = fft2(gb).^kb;
		obj.aux.fgw = fft2(gw).^kw;
		obj.aux.fgh = fft2(gh).^kh;
	case {'krylov'}
		if (obj.bc{1} ~= 'circular' || (obj.bc{2} ~= 'circular'))
			error('only applicable to circular boundary conditions');
		end
		[Dx,Dy,Dxx,Dyy] = setup_matrix();
		%A1 = (I-dt*Dx-dt*Dy);
		%A2 = (I-dt*Dxx-dt*Dyy);
		%L=ichol(A2,struct('type','nofill'));
		A = (I - dt*Dx - dt*Dy - dt*Dxx - dt*Dyy);
		[l,u] = ilu(A,struct('type','nofill'));
	case {'implicit-euler-fourier'}
		[Dx,Dy,Dxx,Dyy] = setup_matrix();
		A = (I - dt*Dx - dt*Dy - dt*Dxx - dt*Dyy);
		x = zeros(obj.n);
		x(1,1)=1;
		x = x(:);
		n = prod(obj.n);
		for idx=1:3
			Ai = A(1+n*(idx-1):n*idx, 1+n*(idx-1):n*idx);
			% impulse response
			r = Ai \ x;
			% Fourier transform
			fr{idx} = fft2(reshape(r,obj.n));
		end
	otherwise
		error('unavailable');
	end % switch diffusion_scheme

	switch (obj.opt.advection_scheme)
	case {'shift'}
		gh  = advection_kernel(obj.dx,dt,obj.n,obj.pmu.vh);
		fgh = fft2(gh);
		obj.aux.fgh = fgh.*obj.aux.fgh;
	case {'advection-diffusion'}
		% nothing to do
	otherwise
		error('not yet implemented');
	end

	y = y0;
	yy(1,:) = y;
	tdxo = 2;
	for tdx=2:nt
		% react half step
		% y = react(t(tdx),y,dt/2);
		y = obj.opt.reaction_scheme(t(tdx),dt/2,y, ...
			@obj.dz_dt_coefficient_react_homogeneous, ...
			@obj.dz_dt_coefficient_react_inhomogeneous);
		%y = obj.opt.reaction_scheme(t(tdx),dt/2,y, @(t,z) obj.dz_dt_coefficient_react_homogeneous(t,z).* z + obj.dz_dt_coefficient_react_inhomogeneous(t,z));

		if (any(y<0))
			% fall back to fdm diffusion, if solution is not smooth
			%y = step_diffuse_fdm_implicit(t,y0,dt,y);
			%warning('solution is not positive');
			error('1 solution is not positive');
			%fallback(tdx) = true;
		end
		% for non-positive reaction, use implicit euler
		%fdx = (y<0);
		%y(fdx) = react_euler_implicit(t(tdx),dt/2,@(t,z) obj.dz_dt_coefficient_react_homogeneous(t,z).* z + obj.dz_dt_coefficient_react_inhomogeneous(t,z));
		% diffuse full step
		switch (obj.opt.diffusion_scheme)
		case {'fdm'}
			y = step_diffuse_fdm_implicit(t(tdx)+dt/2,y,dt);
		case {'fdm-positive','advection-diffusion'}
			y = step_diffuse_fdm_positive(t(tdx)+dt/2,y,dt);
		case {'spectral'}
			y = step_diffuse_analytic(t(tdx)+dt/2,y,dt);
		case {'aid'}
			y = step_diffuse_aid(t(tdx)+dt/2,y,dt);
		case {'krylov'}
			y = step_diffuse_krylov(t(tdx)+dt/2,y,dt);
		case {'implicit-euler-fourier'}
			y = step_diffuse_implicit_euler_fourier(t(tdx)+dt/2,y,dt);
		end % switch diffusion_scheme

		% advect full step
		switch (obj.opt.advection_scheme)
		case {'fmd','aid'}
			% nothing to do, advection is done together with diffusion
		case {'spectral'}
			%y = step_advect_fdm_implicit(t(tdx)+dt/2,y,dt);
			y = step_advect_analytic(t(tdx)+dt/2,y,dt);
		case {'krylov'}
		if (0)
			if (obj.pmu.vh(1) ~= 0)
				y = (I-dt*Dx) \ y;
			end
			if (obj.pmu.vh(2) ~= 0)
				y = (I-dt*Dy) \ y;
			end
		end
		case {'implicit-euler-fourier','advection-diffusion'}
			% nothing to do
		end % switch advection_scheme

		% correct for round off
		y(y<0 & y>-sqrt(eps)) = 0;

		if (any(y<0))
			% fall back to fdm diffusion, if solution is not smooth
			%y = step_diffuse_fdm_implicit(t,y0,dt,y);
			%warning('solution is not positive');
			error('2 solution is not positive');
			%fallback(tdx) = true;
		end
		% react half step
		%y = react(t(tdx)+dt/2,y,dt/2);
		y = obj.opt.reaction_scheme(t(tdx)+dt/2,dt/2,y, ...
			@obj.dz_dt_coefficient_react_homogeneous, ...
			@obj.dz_dt_coefficient_react_inhomogeneous);
		%y = obj.opt.reaction_scheme(t(tdx),dt/2,y, @(t,z) obj.dz_dt_coefficient_react_homogeneous(t,z).* z + obj.dz_dt_coefficient_react_inhomogeneous(t,z));

		% store result
		if (t(tdx) >= to(tdxo))
			yy(tdxo,:) = y;
			tdxo = tdxo+1;
		end % if t(tdx) > to(tdxo)

		if (~isreal(y))
			error('3 solution is not real');
		end
		if (any(y<0))
			% fall back to fdm diffusion, if solution is not smooth
			%y = step_diffuse_fdm_implicit(t,y0,dt,y);
			%warning('solution is not positive');
			error('solution is not positive');
			%fallback(tdx) = true;
		end
	end % for tdx
	yy(end,:) = y;

	function y = step_advect_analytic(t,y,dt)
		c = obj.dz_dt_coefficient(t,y);
		[b,w,h] = obj.extract1(y);
		if (obj.ndim > 1)
			b = reshape(b,obj.n);
			w = reshape(w,obj.n);
			h = reshape(h,obj.n);
		end
		if (0)
			b = advect_analytic(dt,b,obj.L,c{1,2});
			w = advect_analytic(dt,w,obj.L,c{2,2});
		end
			h = advect_analytic(dt,h,obj.L,c{3,2});
	

		if (obj.ndim>1)
			b = flat(b);
			w = flat(w);
			h = flat(h);
		end

		y = [b;w;h];
	end % step_advect_analytic

	function y = step_diffuse_fdm_positive(t,y,dt)
		[b,w,h] = obj.extract1(y);
		if (1 == obj.ndim)
			b = ifft(obj.aux.fgb.*fft(b));
			w = ifft(obj.aux.fgw.*fft(w));
			h = ifft(obj.aux.fgh.*fft(h));
		else
			b = reshape(b,obj.n);
			w = reshape(w,obj.n);
			h = reshape(h,obj.n);
			b = ifft2(obj.aux.fgb.*fft2(b));
			w = ifft2(obj.aux.fgw.*fft2(w));
			h = ifft2(obj.aux.fgh.*fft2(h));
			b = b(:);
			w = w(:);
			h = h(:);
		end
		y = [b;w;h];
	end

	function y = step_diffuse_fdm_implicit(t,y,dt,y0)
		if (nargin()<4)
			y0 = [];
		end
		c  = obj.dz_dt_coefficient(t,y);
		nn = prod(obj.n);
		%[b,w,h] = obj.extract1(y);
		if (1 == obj.ndim)
			error('not yet implemented');
		else
			%D  = c(:,3).*obj.D2;

			A  = (I - (0.5*dt)*D2);
			rhs = y + (0.5*dt)*D2*y;

			if (0)
				y = A\rhs;
			else
			% this is a quick and efficient way to introduce noise
			%y0 = double(single(y));
			maxit = sum(obj.n);
			%icholopt = struct('type','nofill');
			%P = ichol(A,icholopt); 
			P = [];
			[y,flag]  = pcg(A,rhs,pcgtol,maxit,P,P',y0);
			end % else of if 0
			if (flag)
				error('no convergence')
			end
		end % else of 1 == ndim
	end % step_diffuse_fdm_implicit

	function y = step_diffuse_implicit_euler_fourier(t,y,dt,y0)
		[b,w,h] = obj.extract1(y);
		b = reshape(b,obj.n);
		w = reshape(w,obj.n);
		h = reshape(h,obj.n);
		b = ifft2(fr{1}.*fft2(b));
		w = ifft2(fr{2}.*fft2(w));
		h = ifft2(fr{3}.*fft2(h));
		y = [b(:);w(:);h(:)];
	end % step_diffuse_implicit_euler_fourier

	function y = step_diffuse_krylov(t,y,dt,y0)
		if (nargin()<4)
			y0 = [];
		end

		if (1 == obj.ndim)
			% tridiagonal, no iterative solver required
			%y  = (I - dt*Ax) \ y;
			y  = A2 \ y;
		else
			%y = pcg((I-dt*Ax-dt*Ay),y);
			%[y,ncflag] = pcg(A2,y,[],sum(obj.n)); %,L,L');
		%	tic
			[y,ncflag,resn] = cgs(A,y,[],sum(obj.n),l,u);
			%[y,ncflag,resn] = bicg(A,y,[],sum(obj.n),l,u);
			%y = A\y;
			%ncflag = 0;
			%resn   = 0;
			%bbar(end+1,1) = mean(y(1:end/3));
		%	toc
			%tic
			%[y,ncflag] = cgs(A,y,[],sum(obj.n),l,u);
			%toc
			if (ncflag)
				error('pcg did not converge')
			end
		end
	end % step_diffuse_krylov

	function y = step_diffuse_aid(t,y,dt,y0)
		if (nargin()<4)
			y0 = [];
		end

		%[b,w,h] = obj.extract1(y);
		if (1 == obj.ndim)
			%error('not yet implemented');
			y  = (I - dt*Ax) \ y;
		else
			%D  = c(:,3).*obj.D2;
			% note that the matrix can be pre-factored with chol
			if (1 == mod(tdx,2))
				%y  = (I - (0.5*dt)*D2x) \ (y + dt*(0.5*(D2x*y) + 0*D2y*y));
				% impl. euler necessary to assure positivity, trapezoidal is not positive
				% positivity preserving rk? c.f. nusslein 2021, bolley  Crouzeix 1978
				% does the  x,y does matter ? derivatives should commute
				if (prefactor)
					y = ux \ (lx \ y);
					y = uy \ (ly \ y);
				else
					y  = (I - dt*Ax) \ y;
					y  = (I - dt*Ay) \ y;
				end
			else
			%	%y  = (I - (0.5*dt)*D2y) \ (y + dt*(0.5*(D2y*y) + 0*D2x*y));
				if (prefactor)
					y = uy \ (ly \ y);
					y = ux \ (lx \ y);
				else
					y  = (I - dt*Ay) \ y;
					y  = (I - dt*Ax) \ y;
				end
			end
%			tdx
%			min(y)
		end
	end % step_diffuse_aid

	% note that the diffusion coefficient has to be constant in space
	function y = step_diffuse_analytic(t,y,dt)
		c = obj.dz_dt_coefficient(t,y);
		[b,w,h] = obj.extract1(y);

		if (2 == obj.ndim)
			b = reshape(b,obj.n);
			w = reshape(w,obj.n);
			h = reshape(h,obj.n);
		end

		b = diffuse_analytic(dt,b,obj.L,c{1,3});
		w = diffuse_analytic(dt,w,obj.L,c{2,3});
		h = diffuse_analytic(dt,h,obj.L,c{3,3});

		if (2 == obj.ndim)
			b = flat(b);
			w = flat(w);
			h = flat(h);
		end % else of 1 == dim

		y = [b;w;h];
%		y = max(y,0);
	end % step_diffuse_analytic

function [Dx,Dy,Dxx,Dyy] = setup_matrix()
		c  = obj.dz_dt_coefficient(0,y0);
		nn = prod(obj.n);
		c2 = c(:,2);
		c3 = c(:,3);
		I = speye(3*nn);
		Dx = [c2{1}(1)*obj.D1x,            obj.Z, obj.Z;
		                 obj.Z, c2{2}(1)*obj.D1x, obj.Z;
		                 obj.Z,            obj.Z, c2{3}(1)*obj.D1x
	             ];
		Dxx = [c3{1}*obj.D2x,         obj.Z, obj.Z;
		               obj.Z, c3{2}*obj.D2x, obj.Z;
		               obj.Z,         obj.Z, c3{3}(1)*obj.D2x
                     ];
		if (obj.ndim>1)
		Dy = [c2{1}*obj.D1y,        obj.Z, obj.Z;
		              obj.Z,c2{2}*obj.D1y, obj.Z;
		              obj.Z,        obj.Z, c2{3}(2)*obj.D1y
                     ];
		Dyy = [c3{1}*obj.D2y,         obj.Z, obj.Z;
		               obj.Z, c3{2}*obj.D2y, obj.Z;
                               obj.Z,         obj.Z, c3{3}(2)*obj.D2y
                     ];
		else % if ndim > 1
			Dy = [];
			Dyy = [];
		end
end

end % solve_split

