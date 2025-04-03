% 2022-07-05 17:30:32.178128206 +0200
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
%% c.f. koppel
function dz_dt_react = dz_dt_react(obj,t,z)
	% e  : convergence rate
	e = obj.p.e;
	% dm : mortality rate per biomass
	dm = obj.p.dm;
	% km : saturation constant
	km = obj.p.km;
	% Aup : algae density in upper water layer
	Aup = obj.p.Aup;
	% f  : exchange rate between lower and upper layer
	f = obj.p.f;
	% h  : height of lower layer
	h = obj.p.h;
	% c  : consumption rate?
	c = obj.p.c;
	% D  : diffusion rate
	% V  : velocity

	if (isa(Aup,'function_handle'))
		Aup = feval(Aup,t);
	end

	m = z(1:end/2);
	a = z(end/2+1:end);

	dm_dt = e*c*a.*m - dm * km./(km + m).*m;
%	dm_dt = dm_dt + d*(obj.D2x.*m);
	da_dt = (Aup - a)*f - c/h*a.*m;
	da_dt = da_dt + 0.1/sqrt(obj.opt.dt)*randn(size(a));
%	da_dt = da_dt - v*(obj.D1x*a);

	dz_dt_react = [dm_dt; da_dt];
end
 
