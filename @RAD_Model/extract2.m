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
% extract individual variables from joint vector,
% reshape into matrices, if required
%
function varargout = extract2(obj,z)
	varargout = cell(1,obj.nvar);
	[varargout{:}] = obj.extract1(z);

	if (isvector(z))
		if (obj.ndim>1)
			for idx=1:obj.nvar
				varargout{idx} = reshape(varargout{idx},obj.nx);
			end
			%b = reshape(b,obj.n);
			%w = reshape(w,obj.n);
			%h = reshape(h,obj.n);
		end
	else
		nt = size(z,1);
		for idx=1:obj.nvar
			varargout{idx} = reshape(varargout{idx},[nt,rvec(obj.nx)]);
		end
%		b = reshape(b,[nt,obj.n(1),obj.n(2)]);
%		w = reshape(w,[nt,obj.n(1),obj.n(2)]);
%		h = reshape(h,[nt,obj.n(1),obj.n(2)]);
	end
end

