% Mon  5 Jul 17:45:42 CEST 2021
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
%% extract biomass, soil water and surface water from the combined vector
%
function varargout = extract1(obj,z)
	if (ndims(z)>1)
		nt = size(z,1);
	else
		nt = 1;
	end

	n = prod(obj.nx);
	if (isvector(z))
		if (length(z) == obj.nvar)
		for idx=1:obj.nvar
			varargout{idx} = z(idx);
		end
		else	
		for idx=1:obj.nvar
			varargout{idx} = z(1+(idx-1)*n:idx*n);
		end
		end
		%b = z(1:n);	
		%w = z(n+1:2*n);
		%h = z(2*n+1:end);
	else
		for idx=1:obj.nvar
			varargout{idx} = z(:,1+(idx-1)*n:idx*n);
		end
	%	b = z(:,1:n);	
	%	w = z(:,n+1:2*n);
	%	h = z(:,2*n+1:end);
	end
end

