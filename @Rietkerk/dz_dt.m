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
function [dz_dt] = dz_dt(obj,t,z)
	% compute the non-linear flux first, as it is required for the infiltration limiter
	%if (obj.aux.surface_flow)
		if (obj.opt.nonlinear_flow)
		nx     = obj.nx;
		nn     = prod(nx);
		obj.set_drag_coefficient(z);
		dh_dt = obj.aux.zero_inertia.dh_dt(t,z(2*nn+1:end));
		% store rate, as it is co-limiting infiltration
		obj.aux.dh_dt = dh_dt;
		end
		% for nonlinear flow the h-block of AD is zero
		dz_dt = obj.dz_dt_react(t,z) + obj.aux.A*z;
		if (obj.opt.nonlinear_flow)
			dz_dt(2*nn+1:end) = dz_dt(2*nn+1:end) + dh_dt;
		end
	%else % of if surface_flow
	%	dz_dt = obj.dz_dt_react(t,z) + obj.aux.AD_wo_surface_flow*z;
	%end % of else
end

