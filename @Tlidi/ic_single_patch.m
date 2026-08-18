function z = ic_single_patch(obj,sd0)
	L = obj.L;
	[bc,bmax] = obj.homogeneous_state(obj.pmu)
	x = cvec(obj.x);
	z = bc + (bmax-bc)*normpdf(x,L/2+(x(2)-x(1))/2,sd0)/normpdf(0,0,sd0);
end

