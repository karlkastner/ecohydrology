% Wed 29 Jul 15:34:33 CEST 2026
function [z] = random_state(obj)
	dx = obj.dx;
	%p = obj.opt.p_initial;
	p = (1-obj.bmu.mu/obj.pmu.ef);
	class_str = func2str(obj.opt.compute_class);
	switch (obj.ndim)
	case {1}
		z      = obj.pmu.bc*(rand(obj.nx(1),1,class_str)<=p);
		%sqrt(p)*dx(1))/dx(1);
		%z      = obj.pmu.bc*(rand(obj.nx(1),1,class_str)<=sqrt(p)*dx(1))/dx(1);
	case {2}
		z      = obj.pmu.bc*(rand(obj.nx,class_str)<=p);
		% *dx(1)*dx(2))/(dx(1)*dx(2));
		%z      = obj.pmu.bc*(rand(obj.nx,class_str)<=p*dx(1)*dx(2))/(dx(1)*dx(2));
	end
end

