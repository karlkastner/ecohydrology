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
function d2A_dadz = d2A_dadz(obj,t,z,da)
	[b,w,h] = obj.extract1(z);
	nn = prod(obj.nx);
	kb = obj.pmu.kb;
	w0 = obj.pmu.w0;

	d2w_dadb = h./(b + kb) - (h.*(b + kb*w0))./(b + kb).^2;
	d2w_dadh = (b + kb.*w0)./(b + kb);

	d2A_dadz = [sparse([],[],[],nn,3*nn),
		 diag(sparse(d2w_dadb.*da)),sparse([],[],[],nn,nn),diag(sparse(d2w_dadh.*da))
		 diag(sparse(-d2w_dadb.*da)),sparse([],[],[],nn,nn),diag(sparse(-d2w_dadh.*da))
		];
end

