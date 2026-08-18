% Mon 31 May 20:20:46 CEST 2021
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
%% time-derivative of the Rietkerk-pde
%
% function dz_dt = dz_dt(obj,t,z)
% p : parameter vector
% s : standard deviation of paramter
function dz_dt = dz_dt_react(obj,t,z)
	% rate of draining of soil water out into deeper, unmodelled, soil layers
	rdw = obj.drainage_rate(z);

	% rate of soil water evaporation into air
	% does not include transpiration
	re = obj.evaporation_rate(z);

	% rate of soil water uptake (transpiration) by plants
	ru = obj.soil_water_uptake_rate(z);

	% dieback of vegetation
	rdb = obj.dieback_rate(z);

	db_dt = obj.p.cb.*ru - rdb;
	if (obj.aux.isactive(3)) %surface_flow)
		% precipitation rate
		% this must be computed before infiltration, as it limits the infiltrationrate
		% for non-linear infiltration
		rp = obj.precipitation_rate(t,obj.aux.dt);
		obj.aux.rp = rp;

		% rate of infiltration of surface water into soil
		ri = obj.infiltration_rate(z);

		dw_dt = (ri - ru - rdw - re);
		dh_dt = rp - ri;

		% stack output
		dz_dt = [db_dt; dw_dt; dh_dt];
	else
		dw_dt = ( - ru - rdw - re);
		dz_dt = [db_dt; dw_dt];
	end
end % dz_dt_react

