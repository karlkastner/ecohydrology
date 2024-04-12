% Wed 13 Mar 14:55:45 CET 2024
% u_t = u - (1+iA)u*abs(u)^2 + D(u_xx + u_yy)
function dz_dt = dz_dt_react(obj,t,z)
	dz_dt = z - (1.0 + 1i*obj.p.A).*z.*abs(z).^2;
end % dz_dt_react

