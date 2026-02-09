% 2021-07-06 10:56:59.150856622 +0200
% Karl Kästner, Berlin
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
function finish(obj,z)
	% shrink matrices to size containing values
	no                = obj.aux.odx;
	obj.out.to        = obj.out.to(1:no);
	obj.out.zo        = obj.out.zo(1:no,:);
	obj.out.esum      = obj.out.esum(1:no);
	obj.out.n_attempt = obj.out.n_attempt(1:no);
	obj.out.n_error_tolerance_exceeded = obj.out.n_error_tolerance_exceeded(1:no);
	obj.out.n_iter    = obj.out.n_iter(1:no);
	obj.out.n_liter    = obj.out.n_liter(1:no);
	obj.out.n_neg     = obj.out.n_neg(1:no);
	obj.out.n_solver_failed = obj.out.n_solver_failed(1:no);
	obj.out.n_step    = obj.out.n_step(1:no);
	obj.out.runtime   = obj.out.runtime(1:no);
	% TODO z_final and z_last are redundant
	obj.out.z_final = z;
end


