% 2021-07-05 16:48:44.611627340 +0200
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
%% jacobian of the Rietkerk model
%
% dz/dt = A(z)
%
% J = dA(z)/dz
%
function J = jacobian_react(obj,t,z,tflag)
	if (nargin() < 5)
		tflag = false;
	end

	p = obj.p;
	%[b,w,h] = obj.extract1(z);
	if (obj.aux.surface_flow)
		z = reshape(z,[],3);
		h = z(:,3);
	else
		z = reshape(z,[],2);
	end
	b = z(:,1);
	w = z(:,2);

	% derivative of the slowdown term
	if (issym(p.fgb) || p.fgb~=1)
		%s     = (w.^2 + p.fgb*p.kgb.^2)./(w.^2 + p.kgb.^2);
		ds_dw = -(2*p.kgb.^2.*w.*(p.fgb - 1))./(p.kgb.^2 + w.^2).^2;
	else
		ds_dw = 0;
	end

	[U,s]          = obj.soil_water_uptake_rate(z);

	dU_db      = s.*(p.gb*w)./(p.kUw + w);
	dU_dw      = s.*p.gb.*b.*p.kUw./(p.kUw + w).^2;
	dwevap_db  =   p.revap.*p.bevap.*w./(p.bevap+b).^2;
	dwevap_dw  = - p.revap.*p.bevap./(p.bevap+b);
	dwdrain_dw = - p.rdrain.*(w.*(2*p.wdrain + w))./(p.wdrain + w).^2;
	if (obj.aux.surface_flow)
	if (obj.opt.linear_infiltration)
		%dI_db = p.kIb*p.a.*h.*(p.w0 - 1)./(b + p.kIb).^2;
		dI_db = -p.a*(b.^(p.pI-1).*p.kIb^p.pI.*p.pI.*(p.w0 - 1))./((b.^p.pI + p.kIb.^p.pI).^2).*h;
		dI_dh =  p.a*(b.^p.pI + p.kIb.^p.pI*p.w0)./(b.^p.pI + p.kIb.^p.pI);
	else
		if (1 == p.pI)
			dI_db = -p.a*(p.kIb.*(p.w0 - 1))./((b + p.kIb).^2);
		else
			dI_db = -p.a*(b.^(p.pI-1).*p.kIb^p.pI.*p.pI.*(p.w0 - 1))./((b.^p.pI + p.kIb.^p.pI).^2);
		end
		% for constant reate, we have 
		% ri = constant_rate*ie(b)
		% if ri > lim, then there is no dependence on b anymore
		dI_db(obj.aux.ri_limited) = 0;
		dI_dh = zeros(prod(obj.nx),1);
	end
	else
		dI_db = 0;
	end

	if (~issym(z))
	nn = numel(b);
	% TODO use sparse as it is about twice as fast with precomputed indices
	Z = spalloc(nn,nn,0);
	J = [  [diag(sparse(p.cb.*dU_db - s.*p.db)),  diag(sparse(p.cb.*dU_dw + ds_dw.*(p.cb.*U./s - p.db.*b)))];
	       [diag(sparse(- dU_db - dI_db + dwevap_db)),  diag(sparse(-dU_dw - ds_dw.*U./s + dwdrain_dw + dwevap_dw))]];
	if (obj.aux.surface_flow)
	J = [  [J, [Z;
	           diag(sparse(dI_dh))]];
	       [        diag(sparse(dI_db)),           Z, diag(sparse(-dI_dh))]];
	end
	%J = [  [diag(sparse(p.cb.*dU_db - s.*p.db)),  diag(sparse(p.cb.*dU_dw + ds_dw.*(p.cb.*U./s - p.db.*b))),     Z];
	 %      [diag(sparse(- dU_db - dI_db + dwevap_db)),  diag(sparse(-dU_dw - ds_dw.*U./s + dwdrain_dw + dwevap_dw)), diag(sparse(dI_dh))];
	  %     [        diag(sparse(dI_db)),           Z, diag(sparse(-dI_dh))]];
	if (tflag)
		J = J';
	end
	else
	Z = 0;
	J = [   s.*(p.cb*dU_db - p.db),  s.*p.cb.*dU_dw + ds_dw.*(p.cb.*U-p.db.*b),     Z;
	      (- s.*dU_db - dI_db + devap_db),  (-s.*dU_dw - ds_dw.*U + dwdrain_dw + dwevap_dw) , dI_dh;
	               dI_db,           Z, -dI_dh];
	end
end

