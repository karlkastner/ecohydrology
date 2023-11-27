% 2021-07-02 13:29:01.973628215 +0200
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
%% make model parameters symbolic
%
function [p,s] = make_symbolic(obj)
	f_C = fieldnames(obj.pmu);
	p = struct();
	for idx=1:length(f_C)
		p.(f_C{idx}) = sym(f_C{idx});
	end
	obj.p = p;
	obj.aux.D1x = sym('D1x');
	obj.aux.D1y = sym('D1x');
	obj.aux.D2x = sym('D2x');
	obj.aux.D2y = sym('D2y');
	obj.aux.I   = 1;
	obj.aux.Z   = 0;
end

