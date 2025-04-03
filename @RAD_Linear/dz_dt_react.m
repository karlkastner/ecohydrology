% Sun 12 Nov 13:37:24 CET 2023
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
function dz_dt = dz_dt_react(obj,t,z,c)
	nn    = prod(obj.nx);
	z     = reshape(z,nn,obj.nvar);
	if (nargin()<4)
		c= obj.p.c;
	end
	% this is equal to J*z + b, but faster
	dz_dt = (z*c) + obj.p.b';
	dz_dt = reshape(dz_dt,nn*obj.nvar,1);
end

