% 2025-08-05 10:56:47.741732194 +0200
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
function [ru,slowdown] = soil_water_uptake_rate(obj,z)
	p  = obj.p;
	% nb, reshape is faster than extraction
	%z = reshape(z,[],3);
	if (obj.aux.surface_flow)
		z = reshape(z,[],3);
	else
		z = reshape(z,[],2);
	end
	if (p.fgb == 1)
		gb = p.gb;
		slowdown = 1;
	else
		slowdown = (z(:,2).^2 + p.fgb*p.kgb.^2)./(z(:,2).^2 + p.kgb.^2);
		%slowdown = (w.^2 + p.fgb*p.kgb.^2)./(w.^2 + p.kgb.^2);
		gb = slowdown.*p.gb;
	end
	% uptake of water by plants
	ru = gb.*z(:,2)./(z(:,2) + p.kUw).*z(:,1);
end

