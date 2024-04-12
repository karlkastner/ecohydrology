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
function [J] = jacobian_react(obj,t,z,r,outmode,tflag)
	if (nargin()<5)
		outmode = 0;
	end
	if (nargin()<6)
		tflag = false;
	end

	p = obj.p;
	[b,w,h] = obj.extract1(z);

	dU_db = (p.gb*w)./(p.kw + w);
	dU_dw = p.gb*b*p.kw./(p.kw + w).^2;
	dI_db = p.kb*p.a.*h.*(p.w0 - 1)./(b + p.kb).^2;
	dI_dh = p.a.*(b + p.kb*p.w0)./(b + p.kb);


	switch (outmode)
	case {0} % output as matrix
	if (~issym(z))
	nn = prod(obj.nx);
	Z = spalloc(nn,nn,0);
	J = [  [diag(sparse(p.cb*dU_db - p.db)),  diag(sparse(p.cb*dU_dw)),     Z];
	       [diag(sparse(- dU_db - dI_db)),  diag(sparse(-dU_dw - p.rw)), diag(sparse(dI_dh))];
	       [        diag(sparse(dI_db)),           Z, diag(sparse(-dI_dh))]];
	if (tflag)
		J = J';
	end
	else
	Z = 0;
	J = [  p.cb*dU_db - p.db,  p.cb*dU_dw,     Z;
	     - dU_db - dI_db,  -dU_dw - p.rw, dI_dh;
	               dI_db,           Z, -dI_dh];
	end
	case {1}	% output as cell array
		Z = spalloc(obj.nx(1),obj.nx(2),0);
		dU_db = reshape(dU_db,obj.nx);
		dU_dw = reshape(dU_dw,obj.nx);
		dI_db = reshape(dI_db,obj.nx);
		dI_dh = reshape(dI_dh,obj.nx);
		J = { (p.cb*dU_db - p.db), (p.cb*dU_dw), Z; 
	       	      (-dU_db - dI_db),  (-dU_dw - p.rw), dI_dh;
	                 dI_db,           Z, -dI_dh;
                    };	
	case {2} % output as vector-matrix
		% in function jacobian, this is transposed
		if (tflag)
		obj.aux.aa(1,1,:) = (p.cb*dU_db - p.db);
		obj.aux.aa(1,2,:) = (p.cb*dU_dw);
		%obj.aux.aa(1,3,:) = 0; % stays zero
		obj.aux.aa(2,1,:) = (-dU_db - dI_db);
		obj.aux.aa(2,2,:) = (-dU_dw - p.rw);
		obj.aux.aa(2,3,:) = dI_dh;
		obj.aux.aa(3,1,:) = dI_db;
		%obj.aux.aa(3,2,:) = 0;
		obj.aux.aa(3,3,:) = -dI_dh;
		else
		obj.aux.aa(1,1,:) = (p.cb*dU_db - p.db);
		obj.aux.aa(2,1,:) = (p.cb*dU_dw);
		%obj.aux.aa(3,1,:) = 0; % stays zero
		obj.aux.aa(1,2,:) = (-dU_db - dI_db);
		obj.aux.aa(2,2,:) = (-dU_dw - p.rw);
		obj.aux.aa(3,2,:) = dI_dh;
		obj.aux.aa(1,3,:) = dI_db;
		%$obj.aux.aa(2,3,:) = 0; % stays zero
		obj.aux.aa(3,3,:) = -dI_dh;
		end
		J = [];
	end

	if (obj.aux.compute_S)
		r = reshape(r,[],3);
		%d2U_dbdw = p.gb./(p.kw + w) - (p.gb*w)./(p.kw + w).^2;
		%d2U_dw2 = -p.gb*b*p.kw./(p.kw + w).^3;

		obj.aux.S(1,1,:) = r(:,2).*(p.cb.*gb.*kw)./(kw + w)^2;
		obj.aux.S(1,2,:) =  (p.cb.*p.gb.*p.kw.*(p.kw.*r(:,1) - 2*b.*r(:,r2) + r(:,1).*w))./(p.kw + w)^3;
		% S(:,1,3) = 0;
		obj.aux.S(2,1,:) =  (2*a.*h.*p.kb.*r(:,1).*(w0 - 1))./(b + p.kb).^3 - (a.*p.kb.*r(:,3).*(w0 - 1))./(b + p.kb).^2 - (p.gb.*p.kw.*r(:,2))./(p.kw + w).^2;
		obj.aux.S(2,2,:) = -(gb.*kw.*(kw.*r(:,1) - 2.*b.*r(:,2) + r(:,1).*w))./(kw + w).^3;
		obj.aux.S(2,3,:) = -(a.*kb.*r(:,1).*(w0 - 1))./(b + kb)^2;
		obj.aux.S(3,1,:) =  (a.*kb.*(w0 - 1).*(b.*r(:,3) - 2.*h.*r(:,1) + kb.*r(:,3)))./(b + kb).^3;
		% obj.aux.S(:,3,2) = 0;
		obj.aux.S(3,3,:) = (a.*kb.*r(:,1).*(w0 - 1))./(b + kb)^2;
		
	end

%	obj.aux.aa = 0.*obj.aux.aa;
%	J = 0.*J;

if (0)
	if (nargin()<4)
		withda = true;
	end
	U_div_b  = p.gb.*w./(w + p.kw);
        In = p.a.*h.*(b + p.kb.*p.w0)./(b+p.kb);

	if (withda)
		A = [ p.eb*obj.aux.D2x + (p.cb.*U_div_b - p.db).*obj.aux.I, ((b.*p.cb*p.gb*p.kw)./(p.kw + w).^2).*obj.aux.I,     obj.aux.Z
			((p.a.*h)./(b + p.kb) - (p.gb.*w)./(p.kw + w) - (p.a.*h.*(b + p.kb.*p.w0))./(b + p.kb).^2).*obj.aux.I, ...
				p.ew*obj.aux.D2x + (-p.rw  - (b.*p.gb)./(p.kw + w) + (b.*p.gb.*w)/(p.kw + w).^2).*obj.aux.I, ...
				(In./h).*obj.aux.I;
			(p.a.*h.*p.kb.*(p.w0 - 1))./(b + p.kb).^2.*obj.aux.I,    obj.aux.Z, p.eh*obj.aux.D2x + p.vh*obj.aux.D1x - (In./h).*obj.aux.I
		];
	else
		% reaction part only
	A = [ + (p.cb.*U_div_b - p.db).*obj.aux.I, ((b.*p.cb*p.gb*p.kw)./(p.kw + w).^2).*obj.aux.I,     obj.aux.Z
		((p.a.*h)./(b + p.kb) - (p.gb.*w)./(p.kw + w) - (p.a.*h.*(b + p.kb.*p.w0))./(b + p.kb).^2).*obj.aux.I, ...
			(-p.rw  - (b.*p.gb)./(p.kw + w) + (b.*p.gb.*w)/(p.kw + w).^2).*obj.aux.I, ...
			(In./h).*obj.aux.I;
		(p.a.*h.*p.kb.*(p.w0 - 1))./(b + p.kb).^2.*obj.aux.I,    obj.aux.Z, - (In./h).*obj.aux.I
	];
	end
	if (nargout()>1)
	res = [ (b.*p.cb*p.gb.*h.*p.kb)./(h + p.kb).^2
	       b.*( -p.gb.*h.*p.kb./(h + p.kb).^2 + p.a.*w.*p.kw.*(1 - p.w0)./(b + p.kw).^2 )
              -p.R - p.a.*b.*p.kw.*w.*(1 - p.w0)./(b + p.kw).^2];
	end
end
end

