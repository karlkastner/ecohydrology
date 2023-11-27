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
				,'eb', 1*[1,1]     ... % diffusion of biomass (D in klausmeier)
				,'ew', [0,100]     ... % diffusion of water (not in klausmeier)
				,'g',   1          ... % water uptake (R in klausmeier)
				,'l',   1          ... % evaporation (L in clausmeier)
				,'r',   0.077      ... % precipitation (A in Klausmeier)
				,'vw', 182.5*[1,0] ... % water runoff velocity
			);

			obj.pss = struct();
			obj.psl = struct();
			obj.psdist = struct();
			field_C = fieldnames(obj.pmu);
			for idx=1:length(field_C)
				% standard deviation of parameters per unit distance
				% (perturbation once during initialization)
				obj.pss(f{idx}) = 0;
				% spatial correlation lenght of parameters
				obj.psl(f{idx}) = 0;
				% stochastic model (probability distribution) of parameters
				obj.psdist(f{idx}) = [];
			end % for field_C
			
			opt.boundary_condition = {'circular','circular'};
			opt.initial_condition = 'random';

			obj.opt = struct('solver',@ode23 ...
				     ,'dt', 1 ...
				     ,'dto', 100 ...
				     ,'rng', 0 ...
				     , 'path_str', './' ...
				     , 'base_str', 'klausmeier-' ...
				     , 'loadfinal', false ...
				    ); % opt
	
			obj.hashfield_C  = { 'L','dx','T','dt','dto' ...
					,'pmu','pss','psl' ...
					, 'opt.rng' ...
					, 'opt.solver' ...
					, 'initial_condition' ...
					, 'boundary_condition' ...
					}; % hashfield_C
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

