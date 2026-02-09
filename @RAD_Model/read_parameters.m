% 2025-08-02 12:21:02.743509733 +0200
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
function p = read_parameters(obj,filename)
	tab = readtable(filename,'ReadRownames',true,'delimiter',{'\b','\t'});
	p   = rows2vars(tab(:,1));
	% the first column is "originalvariablenames" and is removed here
	p = table2struct(p(:,2:end));
	obj.pmu = copyfields_deep(p,obj.pmu);
end

