% Tue  3 May 09:35:52 CEST 2022
% Karl Kastner, Berlin
%
% rainfall intensity, at which vegetation can only exist through water redistribution
% from bare to vegetated areas, i.e. when it forms patterns
% function Rc = critical_rainfall_depth(obj,p)
function Rc = critical_rainfall_depth(obj,p)
	if (nargin()<2)
		p = obj.p;
	end
	Rc = (p.db*p.kw*p.rw)/(p.cb*p.gb - p.db);
end
