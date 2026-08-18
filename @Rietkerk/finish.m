% 2025-07-23 18:34:37.740271444 +0200
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
function finish(obj,z)
	if (~obj.aux.isactive(3)) %surface_flow)
		z(end+1:end+prod(obj.nx)) = 0;
	end
	finish@RAD_Model(obj,z);
	if (obj.opt.output.store_fluxes)
		no = obj.aux.odx;
		obj.out.dieback = obj.out.dieback(1:no,:);
		obj.out.drainage = obj.out.drainage(1:no,:);
		obj.out.evaporation = obj.out.evaporation(1:no,:);
		obj.out.flow = obj.out.flow(1:no,:);
		obj.out.infiltration = obj.out.infiltration(1:no,:);
		obj.out.precipitation = obj.out.precipitation(1:no,:);
		obj.out.uptake = obj.out.uptake(1:no,:);
		obj.out.zmax = obj.out.zmax(1:no,:);
	end
end

