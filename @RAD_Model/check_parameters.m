% 2025-08-14 11:06:45.839718231 +0200
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
function check_parameters(obj)
	if (       ~isnumeric(obj.T) ...
		|| ~isreal(obj.T) ...
		|| length(obj.T) ~= 2 ...
		|| any(~isfinite(obj.T)))  
		error('T must be a two element real numeric vector with finite values');
	end
%	if (isscalar(obj.T))
%		obj.T = [0,obj.T];
%	end
end
