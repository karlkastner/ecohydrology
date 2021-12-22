% Fri 25 Jun 23:59:47 CEST 2021
function x = solve_stationary(x)
	x = obj.solver_stationary(@obj.stationary_step,x,obj.solver_opt);
end

