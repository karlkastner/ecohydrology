% Mon  2 May 14:18:38 CEST 2022
%
%% evolve the Rietker-PDE in time using a splitting scheme
%
function [t,yy,fallback] = solve_split(obj,to,y0,dt,dmethod)
	AID = 2;
	FDM = 1;
	SPECTRAL = 0;
	if (nargin()<4)
		dmethod = AID;
	end % if nargin < 4
	tol    = 1e-4;
	pcgtol = 1e-5;
	nto = length(to);
	nt  = round((to(end)-to(1))/dt)+1;
	t   = linspace(to(1),to(end),nt)';
	fallback = false(nto,1);
	yy  = zeros(nto,numel(y0));
	prefactor = true;
	resn  =[];
	bbar = [];
	% prepare iteration
	switch (dmethod)
	case(AID)
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
	case (FDM)
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
	case (SPECTRAL)
		% nothing to do, no matrices required
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
%subplot(2,2,1)
%spy(Dyy)
%subplot(2,2,2)
%spy(Dxx)
%full(Dxx(1:obj.n(1),1:obj.n(1)))
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
%rms((Dx(:)))
%rms((Dy(:)))
%rms(flat(Ai-Ai'))
if(0)
subplot(3,3,idx)
imagesc(fftshift(reshape(r,obj.n)))
subplot(3,3,3+idx)
imagesc(fftshift(real(fr{idx})))
colorbar
subplot(3,3,6+idx)
imagesc(fftshift(imag(fr{idx})))
colorbar
end
		end
	otherwise
		error('unavailable');
	end % switch dmethod

	y = y0;
	yy(1,:) = y;
	tdxo = 2;
	for tdx=2:nt
		% react half step
		y = react(t(tdx),y,dt/2);
		% diffuse full step
		switch (dmethod)
		case {FDM}
			y = step_diffuse_fdm_implicit(t(tdx)+dt/2,y,dt);
		case {SPECTRAL}
			y = step_diffuse_analytic(t(tdx)+dt/2,y,dt);
		case {AID}
			y = step_diffuse_aid(t(tdx)+dt/2,y,dt);
		case {'krylov'}
			y = step_diffuse_krylov(t(tdx)+dt/2,y,dt);
		case {'implicit-euler-fourier'}
			y = step_diffuse_implicit_euler_fourier(t(tdx)+dt/2,y,dt);
if (0) %mod(tdx,100)==0)
	imagesc(reshape(y(1:end/3),obj.n));
	title(tdx/nt);
	drawnow
toc
tic
end
		end % switch dmethod

		% advect full step
		switch (dmethod)
		case {FDM,AID}
			% nothing to do, advection is done together with diffusion
		case {SPECTRAL}
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
		case {'implicit-euler-fourier'}
			% nothing to do
		end % switch dmethod
		% react half step
		y = react(t(tdx)+dt/2,y,dt/2);
		% store result
		if (t(tdx) >= to(tdxo))
			yy(tdxo,:) = y;
			tdxo = tdxo+1;
		end % if t(tdx) > to(tdxo)
	end % for tdx
	yy(end,:) = y;

	% TODO this should be better evolved as 3x3 system
	function y = react(t,y,dt)
		% y' = a*y + b
		% y(t) = (exp(a*t)*(b/a + y0) - b/a);
		% fixed point iteration
		% this is better solved with NR, but this requires solution of a linear system
		y0 = y;
		kmax = 100;
		kdx = 0;
		o = ones(prod(obj.n),1);	
		while (1)
		kdx = kdx+1;
		y_ = y;

		% the analytic solution for constant coefficients is
		% y = exp(a*dt)*y0 + (exp(a*dx) - 1)*b/a
		% this solution breaks down for a = 0
		% it is therefore better to approximate the second part as
		% y ~ exp(a*dt)*y0 + dx*b
		% computation time is decreased by approximating the first term as well
		% y ~ exp(a*dt)*y0 + dx*b
		% y ~ (1 + a*dt)*y0 + dx*b
		% the trapezoidal rule does
		% note that the convergence for some reason is not accelerated by aitkens d^2 method
		method = 4;
		switch (method)
		case {1} % implicit euler
		c  = obj.dz_dt_coefficient(t,(y+y0)/2);
		a  = [c{1,1}.*o;
		      c{2,1}.*o;
		      c{3,1}.*o];
		b  = [c{1,4}.*o;
		      c{2,4}.*o;
		      c{3,4}.*o];
%		if (max(abs(a*dt))>1)
%			warning([num2str(max(abs(a*dt))),'max(|a dt|)']);
%		end
		y = (y0 + dt*b)./(1 - 0.5*a*dt);
			
		case {2} % midpoint
		c  = obj.dz_dt_coefficient(t,(y+y0)/2);
		a  = [c{1,1}.*o;
		      c{2,1}.*o;
		      c{3,1}.*o];
		b  = [c{1,4}.*o;
		      c{2,4}.*o;
		      c{3,4}.*o];
		if (max(abs(a*dt))>1)
			warning([num2str(max(abs(a*dt))),'max(|a dt|)']);
		end
		y = ((1 + 0.5*dt*a).*y0 + dt*b)./(1 - 0.5*a*dt);
		%y = exp(dt*a).*y0 + dt*b;
		case {3} % trapezoidal
		c  = obj.dz_dt_coefficient(t,y);
		c0  = obj.dz_dt_coefficient(t,y0);
%		a = c(:,1);
%		b = c(:,4);
		%a0 = c0(:,1);
		%b0 = c0(:,4);
		a  = [c{1,1}.*o;
		      c{2,1}.*o;
		      c{3,1}.*o];
		b  = [c{1,4}.*o;
		      c{2,4}.*o;
		      c{3,4}.*o];
		a0  = [c0{1,1}.*o;
		      c0{2,1}.*o;
		      c0{3,1}.*o];
		b0  = [c0{1,4}.*o;
		      c0{2,4}.*o;
		      c0{3,4}.*o];

		y = ((1 + 0.5*dt*a0).*y0 + 0.5*dt*(b + b0))./(1 - 0.5*a*dt);
		%y = exp(0.5*dt*(a0+a)).*y0 + 0.5*dt*(b + b0);
		case {4}
			% semi-analytic
			c  = obj.dz_dt_coefficient(t,(y0+y)/2);
		%	a  = c(:,1);
		%	b  = c(:,4);
		a  = [c{1,1}.*o;
		      c{2,1}.*o;
		      c{3,1}.*o];
		b  = [c{1,4}.*o;
		      c{2,4}.*o;
		      c{3,4}.*o];
			%y = exp(a*dt)*y0 + (exp(a*dx) - 1)*b/a
			%y = exp(a*dt)*y0 + (exp(a*dx) - 1)*b/a
			y = exp(a*dt).*y0 + dt*b;
			%y = ((1 + 0.5*dt*a).*y0 + dt*b)./(1 - 0.5*a*dt);
		case {5} % euler fw, explicit, no iteration
			c  = obj.dz_dt_coefficient(t,(y0));
			%a  = c(:,1);
			%b  = c(:,4);
			a  = [c{1,1}.*o;
			      c{2,1}.*o;
			      c{3,1}.*o];
			b  = [c{1,4}.*o;
			      c{2,4}.*o;
			      c{3,4}.*o];
			y = ((1 + dt*a).*y0 + dt*b);
		end
		
		dy = rms(y-y_);
		if (dy<tol)
			break;
		end
		if (~isfinite(dy))
			error(['isnan', num2str(kdx), ' ', num2str(dy)]);
		end
		if (kdx>kmax)
			error(['no convergence ', num2str(kdx), ' ', num2str(dy)]);
		end
		end
	end % react

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
			
if (0)
%mod(t+0.5,10) == 0)
subplot(2,2,1)
imagesc(reshape(y(1:end/3),obj.n))
subplot(2,2,2)
plot(resn)
title(t)
drawnow
subplot(2,2,3)
semilogy(bbar)
%pause(0.1)
end
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
		y0 = y;
		c = obj.dz_dt_coefficient(t,y);
		[b,w,h] = obj.extract1(y);

		if (2 == obj.ndim)
			b0 = reshape(b,obj.n);
			w0 = reshape(w,obj.n);
			h = reshape(h,obj.n);
		end

		b = diffuse_analytic(dt,b0,obj.L,c{1,3},false);
		w = diffuse_analytic(dt,w0,obj.L,c{1,3},false);
		h = diffuse_analytic(dt,h,obj.L,c{1,3},false);
%clf
%plot([b0(:,1),b(:,1)])
%hold on
%pause

		if (2 == obj.ndim)
			b = flat(b);
			w = flat(w);
			h = flat(h);
		end % else of 1 == dim

		y = [b;w;h];
		if (any(y<sqrt(eps)))
			% fall back to fdm diffusion, if solution is not smooth
			%y = step_diffuse_fdm_implicit(t,y0,dt,y);
			min(y)
			warning('solution is not positive');
			fallback(tdx) = true;
		end
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

