% Mon 31 May 20:20:46 CEST 2021
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
%% coefficients of the time-derivative of the Rietkerk-pde
%
% function c = dz_dt_coefficient_react_homogeneous(obj,t,z)
function c = dz_dt_coefficient_react_inhomogeneous(obj,t,z)
	if (size(z,2)>1)
		[b,w,h] = obj.extract2(z);
	else
		[b,w,h] = obj.extract1(z);
	end
	if (~isvector(z))
		b=b';
		w=w';
		h=h';
	end

	if (isa(obj.p.R,'function_handle'))
		R = obj.p.R(t,obj.x,obj.y);
	else
		if (isscalar(obj.p.R))
		R = flat(obj.p.R).*ones(prod(obj.n),1);
		else
			R = obj.p.R;
		end
	end

	% infiltration of water into soil
	In = obj.p.a.*obj.infiltration_enhancement(b).*h;

	c = [zeros(prod(obj.n),1);
             In;
	     R;
	     ];
end

