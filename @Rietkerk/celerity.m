% 2021-05-31 21:59:40.374630793 +0200
% Karl Kästner, Berlin
%
%% migration celerity of the pattern
%
% TODO result is awry when pattern is noisy, upwinding necessary ?
%
function [c,cme] = celerity(obj,z)
	if (isvector(z))
		z = rvec(z);
	end
	nt    = size(z,1);
	dz_dt = obj.dz_dt([],z);
	% interestingly, the upwinding seems to be more accurate than
	% central differences, maybe bc upwiding is also used in dz_dt?
	D1x = [obj.aux.D1x,obj.aux.Z,obj.aux.Z;
              obj.aux.Z,obj.aux.D1x,obj.aux.Z;
	      obj.aux.Z,obj.aux.Z,obj.aux.D1x];
	dz_dx = D1x*z';
	if (obj.ndim>1)
	D1y = [obj.aux.D1y,obj.aux.Z,obj.aux.Z;
              obj.aux.Z,obj.aux.D1y,obj.aux.Z;
	      obj.aux.Z,obj.aux.Z,obj.aux.D1y];
	dz_dy = D1y*z';
	end
	
	c   = zeros(nt,obj.ndim);
	cme = zeros(nt,obj.ndim);
	for idx=1:nt
		c(idx,1)   = dz_dx(:,idx) \ dz_dt(:,idx);
		cme(idx,1) = median(dz_dt(:,idx)./dz_dx(:,idx));
		if (obj.ndim>1)
		c(idx,2)   = dz_dy(:,idx) \ dz_dt(:,idx);
		cme(idx,2) = median(dz_dt(:,idx)./dz_dy(:,idx));
		end
	end
end

