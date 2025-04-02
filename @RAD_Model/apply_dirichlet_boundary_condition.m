% set dirichlet BC
function [A,rhs,id] = apply_dirichlet_boundary_condition(obj,A,rhs,val)
	n = obj.nx;
	% TODO 1D
	% TODO non-square domain
	n = n(1);
	% indices of boundary cells
	id = [ (1:n)'
	       (n*(2:n-1))'
	       (n*(3:n))'-1
	       (n-1)*n+(1:n)'];
	nn = n*n;
	% indices for each variable
	% TODO use nvar
	id = [id;
	      id+nn;
	      id+2*nn];
	% indices of matrix diagonals
	ind = sub2ind([3*nn,3*nn],id,id);
	% clear rows
	A(id,:) = 0;
	% set diagonals to 1 for homogeneous bc
	A(ind) = 1;
	% set rhs boundary value
	rhs(id) = val;
end

