% Thu 14 Aug 10:22:57 CEST 2025
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
function z = step_preprocess(obj,t,dt,z)
	% toggle surface flow on and off
	if (obj.aux.surface_flow)
		% check if h is zero
		% TODO no magic numbers	
		hthresh = 1e-4;
		if (   (0 == obj.precipitation_rate(t,dt)) ...
		    && (max(z(2*end/3+1:end))<hthresh) ...
		    )
			obj.aux.surface_flow = 0;
			obj.aux.I = obj.aux.I_without_surface_flow;
			obj.aux.AD = obj.aux.AD_without_surface_flow;
			z = z(1:2*end/3);
		end
	else
		if (obj.precipitation_rate(t,dt) > 0)
			obj.aux.surface_flow = 1;
			z(end+1:(end+prod(obj.nx))) = 0;
			obj.aux.I = obj.aux.I_with_surface_flow;
			obj.aux.AD = obj.aux.AD_with_surface_flow;
		end
	end
end

