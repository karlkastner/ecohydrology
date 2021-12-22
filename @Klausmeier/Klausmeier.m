
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
% c.f. Klausmeier 1999, Regular and Irregular Patterns in Semiarid Vegetation
classdef Klausmeier < handle
	properties
		% runoff velocity
		v

		% rainfall
		r

		% infiltration
		% cancels in non-dimensional model
		
		% die back rate of biomass
		d

		% water-to-biomass conversion rate
		u = 1;

		% first spatial derivative matrix
		D1

		% second spatial derivative matrix
		D2
	end % properties
	methods
		function obj = Klausmeier()
		end % constructor
	end % methods
end % classdef Klausmeier
