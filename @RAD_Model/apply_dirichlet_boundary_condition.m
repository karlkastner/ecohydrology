% 2024-12-13 11:27:38.893271576 +0100
% Karl Kastner, Berlin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%% set dirichlet BC
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

