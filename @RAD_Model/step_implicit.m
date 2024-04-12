% Thu  7 Dec 16:41:17 CET 2023
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
% TODO return iteration info as stat struct, no need for tdx

function [z,stat]= step_implicit(obj,t,z,zold,dt,tt,zz)
	m   = obj.opt.inner_m;
	q   = obj.opt.inner_q;
	stat = struct('rmse',NaN ...
		      ,'t',t ...
		      ,'dt0', NaN ...
		      ,'flag', 0 ...
		      ,'iter', NaN ...
		      ,'ls', struct() ...
		      ,'f0', NaN ...
		      ,'f', NaN  ...
		      ,'inner2', struct() ...
		      ,'flag2', NaN ...
		      ,'fs', NaN ...
		      ,'runtime', NaN ...
		      ,'dt', NaN ...
			);

		switch (obj.opt.inner_solver)
		case {'fminunc'}
			opt = optimoptions('fminunc' ...
					, 'Algorithm','quasi-newton' ...
					, 'SpecifyObjectiveGradient',true ...
					, 'Display','notify' ...
					, 'CheckGradients', false ...
					, 'MaxIterations', maxiter ...
					, 'OptimalityTolerance', obj.opt.inner_tol^2 ...
					, 'FunctionTolerance', obj.opt.inner_tol^2 ...
					, 'HessianApproximation', {'lbfgs', m-1}...
					, 'StepTolerance', 1e-12 ...
					, 'FiniteDifferenceType', 'forward' ...
					);
			%z0 = z;
			[z,stat.fs,stat.exitflag,output,g] = fminunc(@(z1) fun2(t,dt,z1,z), z, opt);
			if (stat.exitflag == 5 || stat.exitflag == 0)
				error(['fminunc ',num2str(stat.exitflag),' ',num2str(stat.fs)])
			end
			obj.out.inner_iter(tdx) = output.iterations;
		case {'quasi-newton'}
			%gamma_ = 1;
			[f_,g_,AtA] = fun2(t,dt,z,z);
			alpha   = max(0,(1.0+sqrt(eps))*(max(sum(abs(AtA),2)./diag(AtA))-2.0));
			options = struct('diagcomp',alpha);
			L  = ichol(AtA,options);
			z0 = z;
			[z,f,stat.exitflag,stat.inner_iter,aw] = lbfgs_man(@(z1) fun2(t,dt,z1,z0), z, m, obj.opt.inner_tol, maxiter, L);
			if (stat.exitflag)
				error('here');
			end
		case {'full-newton'}
			% TODO this requires the d^2r/dz term
			k = 0;
			f    = [];
			inner2 = struc('flag',[],'relres',[],'iter',[]);
			%flag = [];
			timeri = tic();
			while (1)
				[obj.out.f(k+1,1),g,AtA] = fun2(t,dt,z,zold);
				nn_ = length(g);
				if (f(k+1,1) <= tol*tol*nn_)
					break;
				end
				k = k+1;
				switch (obj.opt.innersolver2)
				case {'mldivide'}
					dz = mldivide(AtA,g);
					iter2(k) = 1;
					relres2(k) = 0;
					flag2(k) = 0;
				case {'pcg'}
					alpha   = max(0,(1+sqrt(eps))*(max(sum(abs(AtA),2)./diag(AtA))-2));
					options = struct('diagcomp',alpha);
					L = ichol(AtA,options);
					% do not overwrite f(k) here
					[dz,flag2(k),relres2(k),iter2(k)] = pcg(AtA,g,obj.opt.inner2_tol,obj.opt.inner2_tol,L,L');
					if (obj.out.flag(k,1))
						error('solver did not converge');
					end
				case {'gauss-seidel'}
					x0 = zeros(size(g));
					% do not overwrite f(k) here
					[dz,flag2(k),relres2(k),iter2(k)] = mldivide_gauss_seidel(AtA,g,x0,obj.opt.inner2_tol,opj.opt.inner2_maxiter);
				otherwise
					error('unimplemented method');
				end
				if (~isempty(obj.out.flag) && obj.out.flag(k,1))
					error('inner solver did not convergence');
				end
				dz = -dz;
				% TODO backtracking does not satisfy wolfe condition
				[as,flag(k,2),f(k,2),iter(k,2)] = line_search_backtracking(z,f(k),g, @(z1) fun2(t,dt,z1,zold),dz,1.0);
				z = z + as*dz; 
				if (flag(k,2))
					error('line search did not converge');
				end
			end
			stat.fs             = f(end);
			stat.as = as;
			stat.iter     = size(f,1);
			stat.flag     = flag;
			stat.inner2     = inner2;
		case {'gauss-newton'}
			f    = [];
			k = 0;
			while (1)
				k = k+1;
				stat.iter = k;
				timeri = tic();
				%z    = 0.1*ones(size(z));
				%zold = ones(size(z));
				[stat.f0(k),g,A] = obj.gnfun(t,dt,z,zold);
				% TODO distinguish abstol from reltol
				% TODO convergenc criterion based on gradiend
				if (stat.f0(k) <= obj.opt.inner_tol*obj.opt.inner_tol) % || g'*g < 1e-12)
					break;
				end
				if (k > obj.opt.inner_maxiter)
					warning('innersolver reached maxiter')
					stat.flag(k) = -1;
					return;
				end
				switch (obj.opt.innersolver2)
				case {'bicgstabl'} % newton Krylov
					% preconditioner
					switch (obj.opt.preconditioner)
					case {'ilu'}
						[L,U] = ilu(A);
						prec = {L,U};
					case {'multigrid-java'}
						obj.aux.mg_j.set_coefficients(obj.aux.aa,q*dt*obj.aux.ad,q*dt*obj.aux.e);
						prec = {@mgmfun};
					case {'',[]}
						prec = {};
					otherwise
						error('unimplemented preconditioner')
					end
					%dz = bicgstabl(A,g,obj.opt.inner2_tol,maxiter,prec{:});
					[dz,stat.inner2.flag(k),stat.inner2.relres(k),stat.inner2.iter(k)] = bicgstabl(A,g,obj.opt.inner2_tol,obj.opt.inner2_maxiter,prec{:});
				case {'gmres'}
					[L,U] = ilu(A);
					[dz,stat.inner2.flag(k),stat.inner2.relres(k),stat.inner2.iter(k)] = gmres(A,g,m,obj.opt.inner2_tol,obj.opt.inner2_maxiter,L,U);
					%iter(k,2) = iter_(2);
				case {'gauss-seidel'}
					x0 = zeros(size(g));
					[dz,stat.inner2.flag(k),stat.inner2.relres(k),stat.inner2.iter(k)] = mldivide_gauss_seidel(A,g,x0,obj.opt.inner2_tol,obj.opt.inner2_maxiter);
				case {'mldivide'}
					%dz = (A \ g);
					% spparms('spumoni',2)
					% UMFPACK reorders itselv, colamd makes this slower!
					dz = mldivide(A,g);
					stat.inner2.flag(k) = 0;
					stat.inner2.relres(k) = 0;
					stat.inner2.iter(k) = 1;		
				case {'multigrid'}
					obj.aux.mg_m.set_a(obj.aux.aa);
					obj.aux.mg_m.d = q*dt*rvec(obj.b.ex);
					dz = obj.aux.mg_m.mldivide(g,obj.aux.dz0);
					stat.inner2.flag(k) = obj.aux.mg_m.flag;
					stat.inner2.relres(k)    = obj.aux.mg_m.resn; % TODO relres
					stat.inner2.iter(k) = obj.aux.mg_m.iter;
				case {'multigrid-java'}
%					obj.aux.compute_aa
%					obj.aux.compute_A
					obj.aux.mg_j.set_coefficients(obj.aux.aa,q*dt*obj.aux.ad,q*dt*obj.aux.e);
					gt = reshape(g,obj.nx(1)*obj.nx(2),obj.nvar)';
if (0)
					gt_ = gt;
					g_  = g;
					x = rand(3*prod(obj.nx),1); %size(obj.aux.dz0));
					x(end/3:end)=0;
					xt = reshape(x,obj.nx(1)*obj.nx(2),obj.nvar)';
					res  = A*x - g_;
					res_ = obj.aux.mg_j.resfun_(gt_, xt,0); %obj.aux.dz0);
					res_ = flat(res_');
					%rms(res(end/3:end))
					%rms(res_(end/3:end))
					rms(res(end/3:end)-res_(end/3:end))/rms(res(end/3:end))
					clf
					plot([res])
					hold on
					plot(res_,'--')
					pause
end
					dzt = obj.aux.mg_j.mldivide(gt,obj.aux.dz0);
					dz = flat(dzt');
if (0)
					rms(A*dz - g)
pause 
end
					stat.inner2.flag(k)   = obj.aux.mg_j.flag;
					stat.inner2.relres(k) = obj.aux.mg_j.resn; % TODO relres
					stat.inner2.iter(k)   = obj.aux.mg_j.iter;
				otherwise
					error('unimplemented innersolver2');
				end
				if (stat.inner2.flag(k))
					warning('inner solver 2 did not convergence');
					stat.flag(k) = -2;
					return;
				end
				dz = -dz;
				a_ = 1;
				[as,stat.ls.flag(k),stat.fs(k),stat.ls.iter(k)] = line_search_backtracking(z,stat.f0(k),g, @(z1) obj.gnfun(t,dt,z1,zold),dz,a_);
				z = z + as*dz;
				stat.ls.as(k) = as;
				%stat.fs(k) = f;
				%stat.flag_ls(k) = flag_ls;
				%stat.iter_ls(k) = iter_ls;
				stat.runtime(k) = toc(timeri);
				if (stat.ls.flag(k))
					warning('line search did not converge');
					stat.flag(k) = -3;
					break;
				end
			end % while no convergence
			otherwise
				error('unimplemented solver')
		end % switch opt.innersolver, case newton
	% error estimate
	% TODO length of polynomial and constant according to order of error and scheme
	zz(:,1) = z;
	tt(1)   = t+dt; 
	erel0 = obj.opt.outer_reltol;
	switch (q)
	case {1} % implicit euler
		[e,stat.dt0,eout] = error_step_euler_implicit(tt,zz,erel0,true);
	case {0.5} % trapezoidal
		[e,stat.dt0,eout] = error_step_trapezoidal(tt,zz,erel0,true); 
	otherwise
		% TODO general error estimate for the q-method
		error('q has to be 0.5 or 1 for the error estimate');
	end
	stat.rmse = eout.erms;
	
function x = mgmfun(x)
	x = reshape(x,[],obj.nvar)';
	%x0 = zeros(obj.nvar,prod(obj.nx));
	x = obj.aux.mg_j.cycle1(x);
	%x = obj.aux.mg_j.mldivide(x,x0);
	x=flat(x');
end

% for squared problem
function [f,g,AtA] = fun2(t,dt,z1,z0)
	r = rfun(t,dt,z1,z0);
	% objective
	f = 0.5*(r'*r);
	% gradient
	if (nargout()>1)
		Jt = obj.jacobian(t+dt,z1,r,true);
		A = I - (q*dt)*Jt; 
		g = (A'*r);
		if (nargout()>2)
			AtA = (A'*A);
		end
	end
end % fun2

end % step_implicit

