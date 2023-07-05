% Fri  2 Jul 14:12:49 CEST 2021
% Karl Kästner, Berlin
%
%% homogeneous (not necessarily stable) states of the Rietkerk system
%
function [b,w,h] = homogeneous_state(obj,p,state)
	if (nargin()<2||isempty(p))
		p = obj.pmu;
	end
	if (nargin()<3)
		state = 2;
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
		b = p.cb./p.db.*(p.R.*p.cb.*p.gb - p.db.*(p.R + p.kw.*p.rw)) ...
			         ./(p.cb.*p.gb - p.db);
		w = -(p.db.*p.kw)./(p.db - p.cb.*p.gb).*ones(size(p.R));
		h =  p.R./p.a.*(p.R.*p.cb.*p.db  ...
                           - p.R.*p.cb.^2.*p.gb  ...
                           + p.db.^2.*p.kb  ...
                           - p.cb.*p.db.*p.gb.*p.kb  ...
                           + p.cb.*p.db.*p.kw.*p.rw) ...
		     ./ (p.db.^2.*p.kb.*p.w0 - p.R.*p.cb.^2.*p.gb  ...
                         + p.R.*p.cb.*p.db + p.cb.*p.db.*p.kw.*p.rw  ...
                         - p.cb.*p.db.*p.gb.*p.kb.*p.w0);	
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
end % homogeneous_state

