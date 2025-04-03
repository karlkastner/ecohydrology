% Mon 16 Dec 10:49:12 CET 2024
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
function [c,b,rmse,r2,dz_dt_lin] = statistical_linearization(obj,z,dz_dt_react,ishom)
	nvar = obj.nvar;
	nn = prod(obj.nx);
	z = reshape(z,nn,nvar);
	%C = cov(z)
	d_ = reshape(dz_dt_react,nn,obj.nvar);
	C  = cov(dz_dt_react,d_)

	AA = regression_matrix(z,ishom);

	% get coefficients
	if (1)
		cb = AA \ dz_dt_react;
	else
		W = inv(C);
		% TODO make (nvar*nn)^2
		cb = (obj.aux.A_react*W*obj.aux.A_react) \ (obj.aux.A_react'*W)*dz_dt_react;
	end
	% goodness of fit
	dz_dt_lin = AA*cb;
	rmse = rms(dz_dt_react - dz_dt_lin);
	r2 = 1-rmse.^2./var(dz_dt_react);

	c = reshape(cb(1:nvar^2),obj.nvar,obj.nvar);
	if (~ishom)
	b = cb(obj.nvar^2+1:end);
	else
	b = zeros(obj.nvar,1);
	end

% get covariance matrix
% dz_dt = A (11 .. 1k) (z1) + (b1)
%	    (	...  ) (..) + (..)
%           (k1 .. kk) (zk) + (bk)
function AA = regression_matrix(z,ishom)
	AA = [];
	for idx=1:obj.nvar
		Ai = [zeros(nn,nvar*(idx-1)), z, zeros(nn,nvar*nvar-nvar*(idx))];
		if (~ishom)
			Ai = [Ai,zeros(nn,(idx-1)) ones(nn,1),zeros(nn,nvar-(idx))];
		end
		AA = [AA; Ai];
	end
end

end
