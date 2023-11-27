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
%% time-derivative of the Rietkerk-pde
%
% function dz_dt = dz_dt(obj,t,z)
% p : parameter vector
% s : standard deviation of paramter
% TODO variation of dH is not well implemented good -> define at interfaces between cells
function dz_dt = dz_dt_react(obj,t,z)
	p = obj.p;
%	s = obj.pst;
%%	if (size(z,2)>1)
%		[b,w,h] = obj.extract2(z);
%	else
	[b,w,h] = obj.extract1(z);
%	end
%	if (isvector(z))
%		b = b.';
%		w = w.';
%		h = h.';
%	end
		
	n = prod(obj.nx);

	if (isa(obj.p.R,'function_handle'))
		R = obj.p.R(t);
	else
		R = obj.p.R;
		% isa(obj.p.R,'function_handle'))
	end

	if (isa(obj.p.db,'function_handle'))
		db = obj.p.db(t);
	else
		db = obj.p.db;
	end
	
	% uptake of water by plants
	U = p.gb.*w./(w + p.kw).*b;

	% infiltration of water into soil
	In = p.a.*obj.infiltration_enhancement(b).*h;

	db_dt = p.cb.*U - p.db.*b;
	dw_dt = In - U - p.rw.*w;
	dh_dt = R - In;

	% stack output
	dz_dt = [db_dt; dw_dt; dh_dt];

end % dz_dt_react

