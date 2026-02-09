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
function r = drainage_rate(obj,z)
	p = obj.p;
	%:[b,w,h]=obj.extract2(z);  
	% reshape is faster that extract
	%z = reshape(z,[],3);
	if (obj.aux.surface_flow)
		z = reshape(z,[],3);
	else
		z = reshape(z,[],2);
	end
	%r = p.rdrain.*(w.^2./(p.wdrain + w));
	r = p.rdrain.*(z(:,2).*z(:,2)./(p.wdrain + z(:,2)));
end

