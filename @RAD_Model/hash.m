% Mon  6 Dec 09:43:57 CET 2021
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
%% has the model parameters for filename generation
%
function [key_val, key_str] = hash(obj)
	hashfield_C = obj.hashfield_C;
	% set unused fields to default values
%	f_C = fieldnames(obj.psdist);
%	for idx=1:length(f_C)
%		if (~strcmp(obj.psdist.(f_C{idx}),'geometric-ornstein-uhlenbeck'))
%			obj.psl.(f_C{idx}) = 0;
%		end
%	end
	[key_val, key_str] = hashobj(obj,hashfield_C);
end % hash

