% 2024-12-09 10:04:04.735857954 +0100
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
function sigma = sigma(obj,t,z)
	[b,w,h] = obj.extract1(z);
	if (~isempty(obj.ptsrel))
		sigma = [obj.ptsrel.b.*b;
			 obj.ptsrel.w.*w
			 obj.ptsrel.h.*h
			];
	else
		sigma = [];
	end
	if (~isempty(obj.pts))
		nn = prod(obj.bx);
		sigma(1:nn) = sigma(1:nn) + obj.pts.b;
		sigma(nn+1:2*nn) = sigma(nn+1:2*nn) + obj.pts.w;
		sigma(2*nn+1:end) = sigma(2*nn+1:end) + obj.pts.h;
	end
end

