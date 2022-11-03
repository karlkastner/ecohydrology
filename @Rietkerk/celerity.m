% 2021-05-31 21:59:40.374630793 +0200
% Karl Kästner, Berlin
%
%% migration celerity of the pattern
%
function [c,cme] = celerity(obj,z)
	if (isvector(z))
		z = rvec(z);
	end
	nt    = size(z,1);
	dz_dt = obj.dz_dt([],z);
	% interestingly, the upwinding seems to be more accurate than
	% central differences, maybe bc upwiding is also used in dz_dt?
	D1 = [obj.D1,obj.Z,obj.Z;
              obj.Z,obj.D1,obj.Z;
	      obj.Z,obj.Z,obj.D1];
	dz_dx = D1*z';
	
	c   = zeros(nt,1);
	cme = zeros(nt,1);
	for idx=1:nt
		c(idx)   = dz_dx(:,idx) \ dz_dt(:,idx);
		cme(idx) = median(dz_dt(:,idx)./dz_dx(:,idx));
	end
end

