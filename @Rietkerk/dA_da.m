% 2024-12-12 22:05:38.049672356 +0100
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
function dz_da = dz_da(obj,t,z)
	[b,w,h] = obj.extract1(z);
	% nn = prod(obj.nx);
	nn = numel(b);
	kb = obj.pmu.kb;
	w0 = obj.pmu.w0;
	dw_da = h.*(b + kb.*w0)./(b + kb);
	dz_da = [sparse([],[],[],nn,nn),
		diag(sparse(dw_da)),
		diag(sparse(-dw_da))
		];
end

