% Mon 31 May 20:20:46 CEST 2021
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
function rdb = dieback_rate(obj,z)
	% TODO make mortality a function
	% note that the slow down is for a consistent water balance applied both
	% to the uptate in the biomass and soil water, in contrast to Guttal,
	% who applied it only to the biomass
	if (obj.aux.isactive(3)) %surface_flow)
	%if (obj.aux.surface_flow)
		z = reshape(z,[],3);
	else
		z = reshape(z,[],2);
	end
	if (1 == obj.p.fgb)
		rdb = obj.p.db.*z(:,1);
	else
		slowdown = (z(:,2).^2 + obj.p.fgb*obj.p.kgb.^2)./(z(:,2).^2 + obj.p.kgb.^2);
		rdb = slowdown.*obj.p.db.*z(:,1);
	end
end

