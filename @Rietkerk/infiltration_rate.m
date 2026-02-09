% 2025-08-05 21:27:09.611166686 +0200
% Karl Kastner, Berlin
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
function ri = infiltration_rate(obj,z)
	p = obj.p;
	z = reshape(z,[],3);
	if (obj.opt.linear_infiltration)
		% ri= p.a.*obj.infiltration_enhancement(b).*h;
		ri= p.a.*obj.infiltration_enhancement(z(:,1)).*z(:,3);
	else % constant infiltration
		if (obj.aux.ignore_constant_rates)
			ri = zeros(prod(obj.nx),1);
		else

		% ri = p.a.*obj.infiltration_enhancement(b);
		ri = p.a.*obj.infiltration_enhancement(z(:,1));
		% limit by water availability to preserve positivity
		% no more water can infiltrate than ther is already available (h(0))
		% and being added (removed) by precipitation and surface flow
		h0    = obj.aux.zold(2*end/3+1:end);
		rimax = h0/obj.aux.dt+obj.aux.rp+obj.aux.dh_dt;
		obj.aux.ri_limited = (ri>rimax);
		ri = min(ri,rimax);

		end % else of if ignore_constant rates
	end % else of linear_infiltration
end

