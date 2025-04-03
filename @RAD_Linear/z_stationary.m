% 2024-12-17 17:45:12.951393857 +0100
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
function [z, Ap] = z_stationary(obj,rhs,c,direct)
	nvar = obj.nvar;
	nn = prod(obj.nx);
	A = [];
	cc  = reshape(c(1:nvar^2),nvar,nvar)'; % transpose?
	cd  = c(nvar^2+(1:nvar));
	if (~direct)
		cce = reshape(c(nvar^2+nvar+1:2*nvar^2+nvar),nvar,nvar)';
	end
	D2 = obj.aux.D2x + obj.aux.D2y;
	Ap  = [];
	rhs = reshape(rhs,nn,nvar);
	I   = speye(nn);
	for idx=1:nvar
		Ai = [];
		for jdx=1:nvar
			Aij = cc(idx,jdx)*I;
			if (~direct)
				Aij = diag(sparse(cce(idx,jdx)*e(:,idx)));
			end
			if (idx == jdx)
				Aij = Aij + cd(idx)*D2;
			end
			Ai = [Ai, Aij];
		end
		Ap = [Ap;Ai];
	end
	%if (~direct)
	%	rhs = e(:);
	%end
	%rhs = rhs(:);
	z = Ap \ rhs(:);
end

