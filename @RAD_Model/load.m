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
function [t,y,out] = load(obj)
	out = struct('runtime',[]);
	[oname,oname_final] = obj.filename();
	if (~obj.opt.loadfinal)
		disp(['Loading ',oname]);
		load(oname);
	else
		disp(['Loading ',oname_final]);
		load(oname_final);
		% TODO hot fix
		if (size(t,1) == 2 && size(t,2)>1)
			t = t(1,[1,end])';
		end
		if (size(y,2) == 2)
			y = y';
		end
	end
	obj = copyfields_deep(rad,obj);
end % load

