% 2025-07-28 17:23:55.116572027 +0200
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
function step_postprocess(obj,z)
	% call parent class function
	step_postprocess@RAD_Model(obj,z);

	t = obj.aux.told;
	dt = obj.aux.dt;

	tdx = obj.aux.tdx;	
	odx = obj.aux.odx;

	% TODO this is midpoint integration, implement variantes for trapezoidal
	zmid = 0.5*(z+obj.aux.zold);
	if (obj.opt.nonlinear_flow & obj.aux.surface_flow)
		[b,w,h] = obj.extract1(zmid);
		obj.set_drag_coefficient(b);
		[dh_dt,qi] = obj.aux.zero_inertia.dh_dt(t,h);
		obj.aux.dh_dt = dh_dt;
		if (obj.opt.output.store_step)
			obj.out.step(tdx).boundaryflow(1:2) = [qi(1),qi(end)];
		end
		if (obj.opt.output.store_fluxes)
			obj.out.flow(odx,:) = obj.out.flow(odx,:) + dt*qi;
		end
	end % if obj.opt.nonlinear_flow
	if (obj.opt.output.store_fluxes)
		rp   = obj.precipitation_rate(t,dt);
		obj.aux.rp = rp;
		obj.out.precipitation(odx,1) = obj.out.precipitation(odx,1) + dt*rp;
		rdb = obj.dieback_rate(zmid);
		obj.out.dieback(odx,:) = obj.out.dieback(odx,:) + dt.*rdb';
		u = obj.soil_water_uptake_rate(zmid);
		obj.out.uptake(odx,:) = obj.out.uptake(odx,:) + dt.*u';
		d = obj.drainage_rate(zmid);
		obj.out.drainage(odx,:) = obj.out.drainage(odx,:) + dt.*d';
		e = obj.evaporation_rate(zmid);
		obj.out.evaporation(odx,:) = obj.out.evaporation(odx,:) + dt.*e';
		if (obj.aux.surface_flow)
			ri = obj.infiltration_rate(zmid);
			obj.out.infiltration(odx,:) = obj.out.infiltration(odx,:) + dt.*ri';
		end
		%obj.out.zmax(odx,:)  = max(obj.out.zmax(odx,:),z');
	end % store_fluxes
end % step_postprocess

