% 2025-08-04 13:07:59.239452971 +0200
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
function r = evaporation_rate(obj,z)
	p = obj.p;
	%[b,w,h] = obj.extract2(z);  
	if (obj.aux.surface_flow)
		z = reshape(z,[],3);
	else
		z = reshape(z,[],2);
	end
	%z = reshape(z,[],3);
	% note that reshape only reinterprets but does not copy the data
	%bwh = reshape(z,[],3);
	r = p.revap.*z(:,2).*(p.bevap)./(p.bevap+z(:,1));
	%r = p.revap.*bwh(:,2).*(p.bevap)./(p.bevap+bwh(:,1));
end

