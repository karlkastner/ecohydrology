% Sun 12 Nov 13:33:19 CET 2023
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
%% linear reaction diffusion model
classdef RAD_Linear < RAD_Model
	properties
	end
	methods

		function obj = RAD_Linear(varargin)
			obj = obj.RAD_Linear_();
			if (~isempty(varargin) && isstruct(varargin{1}))
				obj = copyfields_deep(varargin{1},obj);
			else
			    for idx=1:2:length(varargin)-1
				obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			    end
			end
		end
		function obj = RAD_Linear_(obj)
			obj.nvar = 1;
			obj.pmu.ex = 0.1;
			obj.pmu.ey = 0.1;
			obj.pmu.vx = 0;
			obj.pmu.vy = 0;

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

			obj.initial_condition = struct('mu', 1,'sd', 0,'sl',0,'dist',{{'gamma'}});
			obj.boundary_condition = {'circular','circular'};

			obj.opt.rng       = 0;
			obj.opt.solver    = @ode23;
			obj.opt.path_str  = './';
			obj.opt.base_str  = 'rad-linear-';
			obj.opt.loadfinal = false;
			obj.opt.compute_class = @double;
			obj.opt.output_class = 'half';

		end % RAD_Linear_
	end % methods
end % RAD_Linear

