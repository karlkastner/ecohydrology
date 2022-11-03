% 2021-05-31 21:59:40.374630793 +0200
% Karl Kästner, Berlin
%
%% growth rate of biomass of the pattern
%
function r = growth_rate(obj,z)
	if (isvector(z))
		z = rvec(z);
	end
	nt    = size(z,1);
	dz_dt = obj.dz_dt([],z);
	% interestingly, the upwinding seems to be more accurate than
	% central differences, maybe bc upwiding is also used in dz_dt?
	db_dt = dz_dt(1:prod(obj.n),:);
	
	r   = zeros(nt,1);
	for idx=1:nt
		% db_dt = r*b
		r(idx)   = b(:,idx) \ db_dt(:,idx);
		%cme(idx) = median(dz_dt(:,idx)./dz_dx(:,idx));
	end
end

