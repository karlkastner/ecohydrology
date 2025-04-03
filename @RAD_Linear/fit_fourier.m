% 2024-12-17 17:33:25.492004171 +0100
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
% case of direct perturbation of z
% note that when there is an interaction of fz and fe, then this resembles a convolution!!!
function [c,D2f,fz] = fit_fourier(obj,z,e)
	nx = obj.nx;
	nn = prod(nx);
	nvar = obj.nvar;
	z  = reshape(z,[nx,nvar]);
	e  = reshape(e,[nx,nvar]);
	for idx=1:nvar
		fz(:,idx) = flat(fft2(z(:,:,idx)));
		fe(:,idx) = flat(fft2(e(:,:,idx)));
	end

	[fx, fy] = fourier_axis_2d(obj.L,obj.nx);
	ox = 2*pi*fx;
	oy = 2*pi*fy';
	D2f  = -(ox.*ox + oy.*oy);
	D2f = flat(D2f);

	%fz = fz(:);
	% fe1i = [fz1i, .., fzki, 0, ..., 0, -oi^2 fz1i, 0]
	%	 * [a11,..,a1k, .., ak1,..,akk,d1,..,dk]'
	A = [];
	for idx=1:nvar
                Ai = [zeros(nn,nvar*(idx-1)),fz,zeros(nn,nvar*(nvar-idx)), ...
		      zeros(nn,idx-1),D2f.*fz(:,idx),zeros(nn,nvar-idx)];
		A = [A;Ai];
	end
	w = (1:nx)';
%	w = cvec(1./w.^2) + rvec(1./w.^2);
%	w = cvec(1./w.^2) + rvec(1./w.^2);
	w = 1./(w.^2+w'.^2);
	w = repmat(w(:),obj.nvar,1);
	w = diag(sparse(w(:)));
	condest(A'*A)
	% exclude the row for the mean?
	rhs = fe(:);
	%c = A \ rhs;
	c = (A'*w*A) \ (A'*(w*rhs));
	res = A*c - rhs;
	rms(res)/rms(rhs)
end

