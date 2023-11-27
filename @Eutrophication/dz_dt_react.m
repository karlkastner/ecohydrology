% Thu 20 Oct 14:09:22 CEST 2022
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
function dz_dt = dz_dt_react(obj,t,z)
	if (isa(obj.p.a,'function_handle'))
		a = obj.p.a(t);
	else
		a = obj.p.a;
	end
	dz_dt = a - obj.p.b*z + obj.p.r*z.^obj.p.p./(1 + z.^obj.p.q);	
end

