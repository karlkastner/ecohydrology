% Wed 27 Sep 10:15:35 CEST 2023
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
function [t,z,out] = continue_solve(obj,t1,z1,T2)
	z0 = cvec(z1(end,:));
	T_ = [t1(end),T2];
	[t2, z2, out] = obj.solve(z0,T_)
	% concatenating
	t = [cvec(t1);cvec(t2)];
	z = [z1;z2];
end

