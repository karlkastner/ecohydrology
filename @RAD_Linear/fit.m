% Tue 17 Dec 12:33:39 CET 2024
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
function [c,zlin,A] = fit(obj,z,e,rhs,direct)
	nvar = obj.nvar;
	nn = prod(obj.nx);
	z = reshape(z,nn,nvar);
	if (~isempty(e))
	e = reshape(e,nn,nvar);
	end
	D2 = obj.aux.D2x + obj.aux.D2y;
	% note, this fits the coefficients but it does not yield a good estimate
	% build regression matrix
	A = [];
	for idx=1:nvar
		% interaction terms
		% z-terms, D-terms, ze-terms, e0-terms
		Ai = [zeros(nn,nvar*(idx-1)), z, zeros(nn,nvar*(nvar-idx)), ...
		      zeros(nn,(idx-1)), D2*z(:,idx), zeros(nn,nvar-idx)];
		if (~direct)
		      Ai = [Ai, zeros(nn,nvar*(idx-1)), z.*e, zeros(nn,nvar*(nvar-idx)) ];
			%       zeros(nn,idx-1),e,zeros(nn,nvar-idx)];
		end
		A  = [A; Ai];
	end
	if (~direct) % was is direct
		rhs = e;
	end
	% rhs = repmat(e,nvar,1);
	condest(A'*A)
	c   = A \ rhs;
	res = A*c - rhs;
	rms(res)/std(rhs)
	1 - rms(res).^2/var(rhs)

if (1)	% coefficients to z
	% test
	A = [];
	cc  = reshape(c(1:nvar^2),nvar,nvar)'; % transpose?
	cd  = c(nvar^2+(1:nvar));
	if (~direct)
		cce = reshape(c(nvar^2+nvar+1:2*nvar^2+nvar),nvar,nvar)';
	end
	%ce0 = c(2*nvar^2+nvar+1:end);
	Ap = [];
	rhs = reshape(rhs,nn,nvar);
	I = speye(nn);
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
		%rhs(:,idx) = rhs(:,idx) + ce0(idx)*e;
	end
	if (~direct)
		rhs = e(:);
	end
	%rhs = rhs(:);
	zlin = Ap \ rhs(:);
	res  = z(:) - zlin;
	rms(res)/std(z(:))
	1 - rms(res).^2/var(z(:))
end
end

