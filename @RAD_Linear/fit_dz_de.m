% Mon 16 Dec 15:47:54 CET 2024
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
%
% delta_dzdt = dz_dt(e) - dz_dt(0)
function [ce, ce0, rmse, r2] = fit_dz_de(obj, t, z, delta_dzdt, e)
	nn = prod(obj.nx);
	% to determine z :

	% d/de (dz1/dt) = [ c11 e, ..., c1k e] [z1]   [c10 e]
	%              ...                     [..] + [ ... ]
	% d/de (dzk/dt) = [ ck1 e, ..., ckk e] [zk]   [ck0 e]

	% to determine coefficients
	% this is block-diagonal and can be inverted for each variable individually
	% d/de (dzi    = [ ze e, ..., zk e, e] [c11, ..., c1k, c10]' 
	ce = [];
	ce0 = [];
	z = reshape(z,nn,obj.nvar);
	delta_dzdt = reshape(delta_dzdt,nn,obj.nvar);
	for idx=1:obj.nvar
		A   = [z.*e, e];
		cei = (A \ delta_dzdt(:,idx));
		ce  = [ce; rvec(cei(1:end-1))];
		ce0(idx) = cei(end);
		ddzdt_lin(:,idx) = A*cei;
	end
	% goodness of fit
	rmse = rms(ddzdt_lin - ddzdt_lin);
	r2   = 1 - rmse.^2/var(delta_dzdt);
end % fit_dz_de

