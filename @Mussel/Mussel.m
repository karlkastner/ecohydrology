% Sat  9 Nov 10:58:51 CET 2024
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
%
%% c.f. koppel
classdef Mussel < RAD_Model
	properties
	end
	methods
		function obj = Mussel_(obj)
			obj.nvar = 2;
			obj.p = struct();
			% note that l c g eb cancel in the non-dimensional model
			obj.pmu = struct(...
					... % e  : 1, convergence rate
					  'e', 0.2 ...
					... % dm : dm mortality rate per biomass
					, 'dm', 0.02 ...
					... % km : g/m^2 saturation constant
					, 'km', 150 ...
					... % Aup : g/m^3 algae density in upper water layer
					, 'Aup', 1.1 ...
					... % f  : 1/hour, exchange rate between lower and upper layer
					, 'f', 100 ...
					... % h  : m, height of lower layer
					, 'h', 0.1 ...
					... % c m^3/g/hour  : consumption rate?
					, 'c', 0.1 ...
					 ,'ex', [0.005,0] ... % eb,ew,eh, m^2 d^-1
					 ,'ey', [0,0.1] ... % eb,ew,eh, m^2 d^-1
					 ,'vx', [0,360] ... % vb,vw,vh, m d^-1
					 ,'vy', [0,0] ... % vb,vw,vh, m d^-1
			); % struct p

			obj.pss = struct();
			obj.psl = struct();
			obj.psdist = struct();
			field_C = fieldnames(obj.pmu);
			for idx=1:length(field_C)
				% standard deviation of parameters per unit distance
				% (perturbation once during initialization)
				obj.pss.(field_C{idx}) = 0;
				% spatial correlation lenght of parameters
				obj.psl.(field_C{idx}) = 0;
				% stochastic model (probability distribution) of parameters
				obj.psdist.(field_C{idx}) = [];
			end % for field_C
			
			opt.boundary_condition = {'circular','circular'};
			opt.initial_condition = 'random';

			obj.opt.dt = 1; % = struct('solver',@ode23 ...
			obj.opt.dto = 100; %	     ,'dt', 1 ...
			obj.opt.rng = 0; %	     ,'dto', 100 ...
			obj.opt.isreal = true; %	     ,'dto', 100 ...
			%	    		% ,'rng', 0 ...
			%	     , 'path_str', './' ...
			obj.opt.base_str =  'mussel-';
			obj.opt.output_class = @single;
			obj.opt.compute_class = 'double';
			%	     , 'loadfinal', false ...
			%	    ); % opt
	
		end % Mussel_
		function obj = Mussel(varargin)
			obj = obj.Mussel_();
			if (~isempty(varargin) && isstruct(varargin{1}))
				obj = copyfields_deep(varargin{1},obj);
			else
			    for idx=1:2:length(varargin)-1
				obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			    end % for
			end % if
		end % Mussel
	end % properties
end % classdef Mussel

