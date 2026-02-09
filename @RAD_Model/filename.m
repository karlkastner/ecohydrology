% Mon 31 May 20:20:46 CEST 2021
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
function [oname,oname_final,key] = filename(obj)
	key         = obj.hash();
	oname       = [obj.opt.output.path_str,obj.opt.output.base_str,num2str(key),'.mat'];
	oname_final = [obj.opt.output.path_str,obj.opt.output.base_str,num2str(key),'-final.mat'];
end

