% Sat 11 Dec 10:56:37 CET 2021
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
function init(obj)
	% load hashtable
	file_str = [obj.path_str,filesep,obj.map_str];
	mkdir(obj.path_str);
	if (0)
	if (exist(file_str,'file'))
		load(file_str,'map');
		obj.map = map;
	else
		%map = containers.Map('KeyType','double');
		% workaround for Matlab Bug
		map = containers.Map(0,Rietkerk());
		obj.map = map;
	end
	end
end % init

