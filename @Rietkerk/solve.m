% 2021-07-06 10:56:59.150856622 +0200
function [t,zz] = solve(obj, z0)
	T      = linspace(0,obj.T,obj.nt);
	[t,zz] = obj.solver(@obj.dz_dt, T, z0, obj.odeopt);
end

