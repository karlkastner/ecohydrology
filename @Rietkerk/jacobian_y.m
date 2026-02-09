% Thu  7 Dec 09:29:07 CET 2023
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
function J = jacobian_y(obj,t,z,p_react,transpose)
	J = jacobian_y@RAD_Model(obj,t,z,p_react,false);
	if (isfield(obj.opt,'nonlinear_flow') && obj.opt.nonlinear_flow)
		[b,w,h] = obj.extract1(z);
		Jh = obj.aux.zero_inertia.jacobian_j(t,h);
		e3 = sparse(3,3,1,3,3);
		J  = J+kron(e3,Jh);
	end % if ~ linear flow
	if (~transpose)
		J = J';
	end
end
