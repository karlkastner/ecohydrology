% 2022-07-05 17:30:32.178128206 +0200
% c.f. koppel
function dz_dt_mussel(obj,t,z)
	% e  : convergence rate
	% dm : mortality rate per biomass
	% km : saturation constant
	% Aup : algae density in upper water layer
	% f  : exchange rate between lower and upper layer
	% h  : height of lower layer
	% c  : consumption rate?
	% D  : diffusion rate
	% V  : velocity

	m = z(1:end/2);
	a = z(end/2+1:end);

	dm_dt = e*c*a*m - dm * km./(km + m).*m + d*(obj.D2x.*m);
	da_dt = (aup - a)*f - c/h*a*m - v*(obj.D1x*a);
end
 
