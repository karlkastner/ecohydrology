function z = homogenous_state(obj,id)
	p = obj.pmu;
	if (nargin()<2)
		id = 1;

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

	if (isa(Aup,'function_handle'))
		Aup = feval(Aup,0);
	end

	switch (id)
	case {0}
		% trivial equilivbrium
		a = Aup;
		m = 0;
	case {1}
		a = -(dm*km - Aup*e*f*h)/(e*f*h - c*e*km)
 		m = -(f*h*(dm*km - Aup*c*e*km))/(c*(dm*km - Aup*e*f*h))
 		%a = -(dm*km - Aup*e*f*h)/(e*f*h - c*e*km);
		%m = -(dm*f*h*km - Aup*c*e*f*h*km)/(c*dm*km - Aup*c*e*f*h);
 	end
	z = [m;a];

end

