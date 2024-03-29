% Sat 11 Dec 10:56:37 CET 2021
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
%% database for Rietkerk model runs
%
% TODO check for hash collision
% TODO build key table from output files and hashvectors
classdef Rietkerk_Map < handle
	properties
		map
		path_str     = './';
		map_str      = 'rietkerk-map.mat';
		opt          = struct();
	end % properties
	methods
		function obj = Rietkerk_Map(varargin)
			for idx=1:2:length(varargin)
				obj.(varargin{idx}) = varargin{idx+1};
			end
		end
	end % methods
end % classdef

