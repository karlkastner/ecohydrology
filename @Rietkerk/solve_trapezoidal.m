% 2021-07-06 12:45:16.676938153 +0200
% Karl Kästner, Berlin
%
%% trapezoidal time stepping with fixed time step
%
function [t,zz] = solve_trapezoidal(obj,t,z)
	reltol = 1e-6;
	zz = zeros(length(t),length(z));
	z0 = z;
	[A0,i0] = obj.jacobian(t(1),z0);
	I = speye(size(A0));
	dt = t(2)-t(1);
	z1 = z;
	A1 = [];
	i1 = [];
	figure(100)
	clf
	for idx=2:length(t);
		k = 0;
		while (1)
			k = k+1;
			z1_ = z1;
			z1 = step(z1);
			d = rms(z1-z1_)
	%	b = obj.extract2(z1);
	%	figure(100);
	%	subplot(2,3,k)
	%	imagesc(b)
	%	drawnow	
			if (d<reltol)
				break
			end
		end
		[idx k]
		zz(idx,:) = z1;
		z0 = z1;
		A0 = A1;
		i0 = i1;
	%	b = obj.extract2(z1);
	%	figure(100);
	%	subplot(2,3,idx)
	%	imagesc(b)
	%	drawnow	
	end % for idx	

function z1 = step(z1)
	[A1, i1] = obj.jacobian(t(idx),z1);
	z1 = bicgstabl((I - 0.5*dt*A1),z0 + 0.5*dt*A0*z0 + + 0.5*dt*(i0 + i1),[],2e2);% + 0.5*dt*(i0 + i1);
end
end
