% Sun 12 Mar 10:17:12 CET 2023
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
% nb. this is for 1D models
function [z,b,w,h] = initial_condition_periodic(obj,fc,Ric)
	% TODO 2D
	% the simplest way for a consistent condition is to let the foricing
	% term R oscillate and ignore the diffusion and advection terms
	p   = obj.pmu;
	p.db = 0.25;
	p.gb = 0.05;
	p.R = cvec((0.5*Ric).*(1 + cos(2*pi*fc*obj.x)));
	[b,w,h] = obj.homogeneous_state(0,p);
	if (0)
		% TODO 2d
		z = [];
		for idx=1:obj.nvar
			z(:,idx) = a(idx) + s(idx)*sin(2*pi*fc*obj.x) + c(idx)*cos(2*pi*fc*obj.x);
		end
		z = z(:);
	end
	z = [b;w;h];
end

