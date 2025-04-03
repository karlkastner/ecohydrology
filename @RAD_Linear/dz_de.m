% 2024-12-16 17:01:09.145481896 +0100
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
% return the linear part and affine part (residual) of dz_de
function [J, res] = dz_de(obj,t,e)
	ce  = obj.p.ce;
	ce0 = obj.p.ce0;

	% TODO use buffer for construction
	J   = [];
	res = [];
	for idx=1:obj.nvar
	    Ji  = [];
	    res = [res; ce0(idx)*e];
	    for jdx=1:obj.nvar
		Ji = [Ji, diag(sparse(ce(idx,jdx)*e))];
	    end
	    J = [J; Ji];
	end
end % dz_de

