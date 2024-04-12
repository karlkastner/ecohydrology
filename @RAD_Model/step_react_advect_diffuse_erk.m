% Wed 13 Mar 13:34:35 CET 2024
function [z, stat] = step_react_advect_diffuse_erk(obj,t,z,dt,varargin)
	dx   = obj.dx;
	a    = [obj.pmu.vy; obj.pmu.vx];
	e    = [obj.pmu.ey; obj.pmu.ex];
	nfun = @(z) obj.dz_dt_react(t,z);

	z = step_react_advect_diffuse_erk(dt,dx,a,e,nfun,reshape(z,[obj.nx,obj.nvar]),obj.opt.isreal);
	% TODO
	stat = struct('rmse',NaN,'dt0',dt,'flag',0);
end

