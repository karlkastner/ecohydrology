function J = jacobian_react(obj,z)
	p = obj.pmu;

	% e  : convergence rate
	e = p.e;
	% dm : mortality rate per biomass
	dm = p.dm;
	% km : saturation constant
	km = p.km;
	% Aup : algae density in upper water layer
	Aup = p.Aup;
	% f  : exchange rate between lower and upper layer
	f = p.f;
	% h  : height of lower layer
	h = p.h;
	% c  : consumption rate?
	c = p.c;
	% D  : diffusion rate
	% V  : velocity

	m = z(1);
	a = z(2);
	J = [a*c*e - (dm*km)/(km + m) + (dm*km*m)/(km + m)^2,         c*e*m
                                       -(a*c)/h, - f - (c*m)/h];
end
