% 2023-10-18 12:55:55.325396302 +0200
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
%
% flag is true when test fails
function [flag,z0,dz_dt] = test_homogeneous_state()
	rad = Dodorico();
	z0 = rad.homogeneous_state();
	rad.p = rad.pmu;
	dz_dt = rad.dz_dt_react(0,z0);
	flag = max(abs(dz_dt)>sqrt(eps))
end

