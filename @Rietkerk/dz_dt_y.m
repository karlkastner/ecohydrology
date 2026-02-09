% Wed 18 Oct 08:42:39 CEST 2023
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
function [dz_dt] = dz_dt_x(obj,t,z,p,varargin)
	dz_dt = (p*obj.dz_dt_react(t,z,varargin{:}) + obj.aux.Ay*z);
	if (isfield(obj.opt,'nonlinear_flow') && obj.opt.nonlinear_flow)
		nx     = obj.nx;
		nn     = prod(nx);
		[b,w,h]=obj.extract1(z);
		dh_dt = obj.aux.zero_inertia.dh_dt_j(t,h);
		dz_dt(2*nn+1:end) = dz_dt(2*nn+1:end) + dh_dt;
	end % if ~ linear flow
end

