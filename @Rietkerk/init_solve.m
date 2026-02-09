% 2025-08-06 19:13:15.897720793 +0200
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
function init_solve(obj)
	init_solve@RAD_Model(obj);
	if (isfield(obj.aux,'AD'))
		obj.aux.AD_with_surface_flow = obj.aux.AD;
		obj.aux.I_with_surface_flow  = obj.aux.I;
		obj.aux.AD_without_surface_flow = obj.aux.AD(1:2*end/3,1:2*end/3);
		obj.aux.I_without_surface_flow  = obj.aux.I(1:2*end/3,1:2*end/3);
	end

	obj.aux.ignore_constant_rates = false;
	obj.aux.surface_flow = true;
	if (obj.opt.output.store_fluxes)
		no = length(obj.out.to);
		output_class = func2str(obj.opt.output_class);
		obj.out.dieback      = zeros(no,prod(obj.nx),output_class);
		obj.out.drainage     = zeros(no,prod(obj.nx),output_class);
		obj.out.evaporation  = zeros(no,prod(obj.nx),output_class);
		% TODO overflow for half
		obj.out.flow         = zeros(no,prod(obj.nx+1),'single');
		obj.out.infiltration = zeros(no,prod(obj.nx),output_class);
		obj.out.precipitation = zeros(no,1,output_class);
		obj.out.uptake       = zeros(no,prod(obj.nx),output_class);
		obj.out.zmax         = zeros(no,3*prod(obj.nx),output_class);
	end	
	if (isstruct(obj.pmu.R))
		R = obj.pmu.R;
		t = obj.pmu.R.time;
		t = t-t(1);
		depth = obj.pmu.R.depth;

		[Ri,tevent] = rate_interpolator_from_time_series(t,depth);

		obj.p.R = Ri;

		% find precipitation event times
		switch (obj.opt.event)
		case {'none'}
			obj.aux.tevent = [0,inf];
		case {'start'}
			tflag = ((depth>0) & [true; 0==depth(1:end-1)]);
			obj.aux.tevent = [t(tflag);inf];
		case{'change'}
			tflag = ([true; depth(2:end) ~= depth(1:end-1)]);
			obj.aux.tevent =[ t(tflag); inf];
		otherwise
			error('here');
		end % switch obj.opt.event
	end % isstruct obj.pmu.R
end

