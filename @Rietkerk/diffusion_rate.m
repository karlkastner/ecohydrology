% 2022-06-26 16:39:01.630950652 +0200
% Karl Kästner, Berlin
%
% estimate the diffusion rate e of the pattern
%
% under the assumption that diffusion is independent of reaction and advection:
%
% dz/dt = e*(d^2/dx^2 + d^2/dy^2) z

function e = diffusion_rate(obj,z)
	if (isvector(z))
		z = rvec(z);
	end
	nt    = size(z,1);
	dz_dt = obj.dz_dt([],z);

	D2 = [obj.D2,obj.Z,obj.Z;
              obj.Z,obj.D2,obj.Z;
	      obj.Z,obj.Z,obj.D2];
	d2z_dx2 = D2*z';
	
	e   = zeros(nt,1);
%	cme = zeros(nt,1);
	for idx=1:nt
		e(idx)   = d2z_dx2(:,idx) \ dz_dt(:,idx);
%		cme(idx) = median(dz_dt(:,idx)./dz_dx(:,idx));
	end
end % diffusion_rate

