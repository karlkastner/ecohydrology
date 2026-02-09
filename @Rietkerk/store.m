% 2025-07-23 18:34:37.740271444 +0200
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
function store(obj,t,z)
	% call superclass storage function
	store@RAD_Model(obj,t,z);

	if (obj.opt.output.store_fluxes)
		% reallocate
		no = size(obj.out.uptake,2);
		if (obj.aux.odx+1>no)
			no = 2*no;
			% expand arrays
			obj.out.dieback(no,1) = 0;
			obj.out.drainage(no,1) = 0;
			obj.out.evaporation(no,1) = 0;
			obj.out.flow(no,1) = 0;
			obj.out.infiltration(no,1) = 0;
			obj.out.precipitation(no,1) = 0;
			obj.out.uptake(no,1) = 0;
			obj.out.zmax(no,1) = 0;
		end
	end
end

