% Fri  2 Jul 14:12:49 CEST 2021
%
% homogeneous (not necessarily stable) states of the rietkerk system
%
function [b,w,h] = homogeneous_state(obj,p,state)
	if (nargin()<2||isempty(p))
		p = obj.pmu;
	end
	if (nargin()<3)
		state = 1;
	end
	switch(state)
	case{0}
		% unvegetated
		b = 0;
       	 	w = p.R./p.rw;
		h = p.R./(p.a.*p.w0);
	case{1}
		% state 1, vegetated (when R sufficiently large),
		%          otherwise biomass is negative
		b = p.cb./p.db.*(p.R.*p.cb.*p.gb - p.db.*(p.R + p.kw.*p.rw)) ...
			         ./(p.cb.*p.gb - p.db);
		w = -(p.db.*p.kw)./(p.db - p.cb.*p.gb);
		h =  p.R./p.a.*(p.R.*p.cb.*p.db  ...
                           - p.R.*p.cb.^2.*p.gb  ...
                           + p.db.^2.*p.kb  ...
                           - p.cb.*p.db.*p.gb.*p.kb  ...
                           + p.cb.*p.db.*p.kw.*p.rw) ...
		     ./ (p.db.^2.*p.kb.*p.w0 - p.R.*p.cb.^2.*p.gb  ...
                         + p.R.*p.cb.*p.db + p.cb.*p.db.*p.kw.*p.rw  ...
                         - p.cb.*p.db.*p.gb.*p.kb.*p.w0);	

	end % state
end % homogeneous_state

