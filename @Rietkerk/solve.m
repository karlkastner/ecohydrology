% 2021-07-06 10:56:59.150856622 +0200
% Karl Kästner, Berlin
%
%% solve the Rietkerk model
%
function [t,zz] = solve(obj, z0, T)
	opt = obj.odeopt;
	if (nargin()<2)
		z0 = obj.z0;
	end
	if (nargin()<3)
		nto    = max(2,round(obj.T/obj.dto));
		T      = linspace(0,obj.T,nto)';
	end
	switch (func2str(obj.opt.solver))
	%case {'ode23t','ode23tb'}
	case {'solve_split','Rietkerk.solve_split'}
		%opt = odeset(opt,'Jacobian',@obj.jacobian);
		% TODO no magic numbers
		%[t,zz] = obj.solve_split(@obj.dz_dt, T, z0, method); %, obj.odeopt);
		[t,zz] = obj.solve_split(T, cvec(z0));
	otherwise
		[t,zz] = obj.opt.solver(@obj.dz_dt, T, z0, obj.odeopt);
	end
end

