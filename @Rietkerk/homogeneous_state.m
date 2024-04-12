% Fri  2 Jul 14:12:49 CEST 2021
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
%% homogeneous (not necessarily stable) states of the Rietkerk system
%
function [b,w,h,J,v,e] = homogeneous_state(obj,p,state,outz)
	if (nargin()<2||isempty(p))
		p = obj.pmu;
	end
	if (nargin()<3)
		state = 2;
	end
	if (isnumeric(p.db))
		db = p.db;
	else
		db = p.db(0);
	end
	switch(state)
	case{0}
		% unvegetated
       	 	w = p.R./p.rw;
		h = p.R./(p.a.*p.w0);
		b = zeros(size(w));
	case{1}
		% state 1, vegetated (when R sufficiently large),
		%          otherwise biomass is negative
		b = p.cb./db.*(p.R.*p.cb.*p.gb - db.*(p.R + p.kw.*p.rw)) ...
			         ./(p.cb.*p.gb - db);
		w = -(db.*p.kw)./(db - p.cb.*p.gb).*ones(size(p.R));
		h =  p.R./p.a.*(p.R.*p.cb.*db  ...
                           - p.R.*p.cb.^2.*p.gb  ...
                           + db.^2.*p.kb  ...
                           - p.cb.*db.*p.gb.*p.kb  ...
                           + p.cb.*db.*p.kw.*p.rw) ...
		     ./ (db.^2.*p.kb.*p.w0 - p.R.*p.cb.^2.*p.gb  ...
                         + p.R.*p.cb.*db + p.cb.*db.*p.kw.*p.rw  ...
                         - p.cb.*db.*p.gb.*p.kb.*p.w0);	
	case {2}
		% state dependent on (local) water availability
		Rc         = obj.critical_rainfall_depth(p);
		[b,w,h]    = obj.homogeneous_state(p,0);
		[b1,w1,h1] = obj.homogeneous_state(p,1);
		fdx     = p.R > Rc;
		b(fdx)  = b1(fdx);
		w(fdx)  = w1(fdx);
		h(fdx)  = h1(fdx);
	otherwise
		error('option unavailable');
	end % state
	if (nargout()>3)
		J = obj.jacobian(0,[b,w,h],false);
		[v,e] = eig(J);
	end
	if (nargin()>3&&outz)
		o = ones(prod(obj.nx),1);
		b = [b.*o;w.*o;h.*o];
	end
end % homogeneous_state

