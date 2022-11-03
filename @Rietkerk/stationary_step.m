% 2021-07-05 16:49:05.311935922 +0200
% Karl Kästner, Berlin
%
%% quasi-stationary time-step
%
function z = stationary_step(obj,z)
	[A,rhs]  = obj.jacobian(obj,t,z);
	%z  = minres(A,z);
	z = bicgstabl(A,rhs,[],maxit,[],[],x);
end

