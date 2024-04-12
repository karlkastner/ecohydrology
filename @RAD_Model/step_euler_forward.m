% Mon 16 Oct 09:40:37 CEST 2023
function [z,stat] = step_euler_forward(obj,t,z,dt,tt,zz)
	dz_dt = obj.dz_dt(t,z);
	z = z + dt*dz_dt;
	stat = struct('rmse',NaN,'dt0',NaN,'flag',0);
end

