% 2021-10-27 19:38:39.504517593 +0200
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
%% infiltration enhancement of the Rietkerk model
%
function ie = infiltration_enhancement(obj,b)
	%ie = (b + obj.p.kIb.*obj.p.w0)./(b+obj.p.kIb);
	if (1 == obj.p.pI)
		ie = (b + obj.p.kIb*obj.p.w0)./(b + obj.p.kIb);
	else
		ie = (b.^obj.p.pI + obj.p.kIb.^obj.p.pI*obj.p.w0)./(b.^obj.p.pI + obj.p.kIb.^obj.p.pI);
	end
end


