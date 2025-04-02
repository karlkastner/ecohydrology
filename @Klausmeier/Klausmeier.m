% 2021-12-22 17:03:30.000000000 +0100
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
%
% c.f. Klausmeier 1999, Regular and Irregular Patterns in Semiarid Vegetation
classdef Klausmeier < RAD_Model
	properties
		noise
	end
	methods
		function obj = Klausmeier_(obj)
			obj.nvar = 2;
			obj.p = struct();
			% note that l c g eb cancel in the non-dimensional model
			obj.pmu = struct(... 
				 'c',   1          ... % conversion rate of water to biomass (J in Klausmeier)
				,'d',  0.045       ... % death rate of biomass (M in klausmeier)	
				,'g',   1          ... % water uptake (R in klausmeier)
				,'l',   1          ... % evaporation (L in clausmeier)
				,'r',   0.077      ... % precipitation (A in Klausmeier)
				,'ex', [1,0]	   ... % diffusion of biomass (D) and water (0 in klausmeier)
				,'ey', [1,0]	   ... %
				,'vx', [0,182.5]   ... % advection of biomass (0) and water
				,'vy', [0,0]   ... % advection of biomass (0) and water
			);

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
				obj.pst.(field_C{idx}) = 0;
				% stochastic model (probability distribution) of parameters
				obj.psdist.(field_C{idx}) = [];
			end % for field_C
			
			obj.opt.boundary_condition = {'circular','circular'};
			obj.opt.initial_condition = 'random';

			%obj.obj.opt = struct('solver',@ode23 ...
			obj.opt.dt = '1';
		        obj.opt.dto = 100;
		        obj.opt.rng = 0;
		        obj.opt.path_str = './';
		        obj.opt.base_str = 'klausmeier-';
		        obj.opt.loadfinal = false;
		        obj.opt.isreal = true;
	
		end % klausmeier_
		function obj = Klausmeier(varargin)
			obj = obj.Klausmeier_();
			if (~isempty(varargin) && isstruct(varargin{1}))
				obj = copyfields_deep(varargin{1},obj);
			else
			    for idx=1:2:length(varargin)-1
				obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			    end
			end
		end % constructor
%		function [x,y] = x(obj)
%			dx = obj.L./obj.nx;
%			x = (0:obj.nx(1)-1)*obj.L(2)/obj.nx(1);
%			if (nargout()>1)
%			y = (0:obj.nx(2)-1)*obj.L(2)/obj.nx(2);
%			end
%		end
%		function dx = dx(obj)
%			dx = obj.L./obj.dx;
%		end
	end % methods
end % classdef Klausmeier

