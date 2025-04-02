% 2021-05-31 21:59:40.374630793 +0200
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
%% migration celerity of the pattern
%
% TODO result is awry when pattern is noisy, upwinding necessary ?
%
function [c,cme] = celerity(obj,z,individual)
	if (nargin()<3)
		individual = 0;
	end
	if (~isvector(z))
		nt    = size(z,1);
	else
		nt = 1;
	end
	if (isvector(z))
		z = cvec(z);
	end
	if (~isfield(obj.aux,'AD'))
		obj.init_advection_diffusion_matrix();
	end
	z = double(z);
	dz_dt = obj.dz_dt([],z);
	% interestingly, the upwinding seems to be more accurate than
	% central differences, maybe bc upwiding is also used in dz_dt?
	n = prod(obj.nx);
	obj.aux.Z = sparse(n,n);
	D1xl = [obj.aux.D1xl,obj.aux.Z,obj.aux.Z;
              obj.aux.Z,obj.aux.D1xl,obj.aux.Z;
	      obj.aux.Z,obj.aux.Z,obj.aux.D1xl];
	dz_dx = D1xl*z;
	if (obj.ndim>1)
	D1yl = [obj.aux.D1yl,obj.aux.Z,obj.aux.Z;
              obj.aux.Z,obj.aux.D1yl,obj.aux.Z;
	      obj.aux.Z,obj.aux.Z,obj.aux.D1yl];
	dz_dy = D1yl*z;
	end
	
	c   = zeros(nt,obj.ndim);
	cme = zeros(nt,obj.ndim);
	if (~individual)
	for idx=1:nt
		cme(idx,1) = median(dz_dt(:,idx)./dz_dx(:,idx));
		if (obj.ndim==1)
		c(idx,1)   = dz_dx(:,idx) \ dz_dt(:,idx);
		else
		c(idx,1:2)   = [dz_dx(:,idx),dz_dy(:,idx)] \ dz_dt(:,idx);
		cme(idx,2) = median(dz_dt(:,idx)./dz_dy(:,idx));
		end
	end
	else
	nn = prod(obj.nx);
	for idx=1:nt
		for jdx=1:obj.nvar
		if (obj.ndim==1)
		c(idx,1,jdx)   = dz_dx(1+(jdx-1)*nn:jdx*nn,idx) \ dz_dt(1+(jdx-1)*nn:jdx*nn,idx);
		else
		c(idx,1:2,jdx)   = [dz_dx(1+(jdx-1)*nn:jdx*nn,idx) dz_dy(1+(jdx-1)*nn:jdx*nn,idx)] \ dz_dt(1+(jdx-1)*nn:jdx*nn,idx);
		end
		%cme(idx,1) = median(dz_dt(:,idx)./dz_dx(:,idx));
		%c(idx,2)   = dz_dy(:,idx) \ dz_dt(:,idx);
		%cme(idx,2) = median(dz_dt(:,idx)./dz_dy(:,idx));
		end
	end
	end
end

