% Tue  3 May 09:35:52 CEST 2022
% Karl Kastner, Berlin
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
% rainfall intensity, at which vegetation can only exist through water redistribution
% from bare to vegetated areas, i.e. when it forms patterns
% function Rc = critical_rainfall_depth(obj,p)
function Rc = critical_rainfall_depth(obj,p)
	if (nargin()<2)
		p = obj.p;
	end
	Rc = (p.db*p.kUw*p.revap)/(p.cb*p.gb - p.db);
	%Rc = (p.db*p.kw*p.rw)/(p.cb*p.gb - p.db);
end
