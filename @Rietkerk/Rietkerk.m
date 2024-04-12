% Mon  5 Jul 18:00:15 CEST 2021
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
%% c.f. Rietkerk et al. 2002, Self-Organization of Vegetation in Arid Ecosystems 
%
classdef Rietkerk < RAD_Model
	properties
		% all inherited from RAD_Model
	end % properties
	methods
	% pseudo-constructor
	function obj = Rietkerk_(obj)
		obj.nvar = 3;
		% mean of the physical parameters
		obj.pmu = struct(         ...
			  'cb',	10    ... % g mm^-1 m^-2
			... % range given as 0 to 0.5 by rietkerk
			 ,'db',	 0.25 ... % d^-1
			 ,'kw',	 5    ... % mm
			 ,'gb',	 0.05 ... % mm g^-1 m^-2 d^-1
			 ,'a',	 0.2  ... % d^-1
			 ,'rw',	 0.2  ... % d^-1
			 ,'ex', [0.1,0.1,100] ... % eb,ew,eh, m^2 d^-1
			 ,'ey', [0.1,0.1,100] ... % eb,ew,eh, m^2 d^-1
			 ,'vx', [0,0,0] ... % vb,vw,vh, m d^-1
			 ,'vy', [0,0,0] ... % vb,vw,vh, m d^-1
			 ...,'eb',	 0.1  ... % m^2 d^-1
			 ...,'ew',	 0.1  ... % m^2 d^-1
		 	 ...,'eh',	 100    ... % m^2 d^-1
			 ...,'vh',	 0.0  ... %  
			 ...% range given as 0 to 3 by Rietkerk
			 ,'R',	 1    ... % mm d^-1
			 ,'w0',	 0.2  ...   % 1
			 ,'kb',	 5 ...     % g m^-2
			 ... , 'Manning', 0.055 ... % c.f. caviedes 2022 
		         ... , 'dzb_dx', 0 ...
		); % struct p

		field_C = fieldnames(obj.pmu);
		for idx=1:length(field_C)
			% standard deviation of parameters per unit distance
			% (perturbation once during initialization)
			obj.pss.(field_C{idx}) = 0;
			% spatial correlation lenght of parameters
			obj.psl.(field_C{idx}) = 0;
			% stochastic model (probability distribution) of parameters
			obj.psdist.(field_C{idx}) = [];
		end % for 

		obj.initial_condition = struct('mu',[0,0,0] ...
					       ,'sd',[0,0,0] ...
					       ,'sl',[0,0,0]);
		obj.initial_condition.dist = {[],[],[]}; 

		% boundary conditions
		obj.boundary_condition = {'circular','circular'};

		% derivative matrices
		obj.aux      = struct('fgb',1,'fgw',1,'fgh',1);

		obj.hashfield_C  = { 'L','dx','T','opt.dt','opt.dto' ...
				,'pmu','pss','psl' ...
				,'initial_condition' ...
				,'boundary_condition' ...
				,'opt.rng' ...
				,'opt.solver' ...
				... ,'opt.reaction_scheme' ...
				... ,'opt.advection_scheme' ...
				... ,'opt.diffusion_scheme' ...
				,'opt.compute_class' ...
				}; % hashfiled_C

		%obj.opt = struct(  ...
                             ...  'reaction_scheme',  @step_react_ralston_exp ...
			     ... , 'diffusion_scheme', 'advection-diffusion' ...
			     ... , 'advection_scheme', 'advection-diffusion' ...
			     ... , 'tol_implicit_iteration', 1e-4 ...
			     ... , 'tol_pcg', 1e-5 ...
			     ... , 'solver_stationary', @picard ...
		%obj.opt.solver = 'solve_split' ...
		%obj.opt.rng = 0; ... % state of random number generator at start
		%obj.opt.zero_inertia', 'false' ...
			     .... , 'order', 2 ...
			     ... , 'symbolic', false ...
			     ... , 'heterogeneity_model', 'geometric-ornstein-uhlenbeck' ...
		%obj.opt.heterogeneity_integration_order = 9; % ...
		%	     , 'compute_class', @single ...
		%	     , 'output_class', 'half' ...
		%	     , 'adaptive_time_step', false ...
			     ... , 'dytol', 1e-2 ...
		%	     , 'path_str', './' ...
		obj.opt.base_str = 'rietkerk-';
		obj.opt.isreal = true;
		%loadfinal', false ...
		             % opt
		% solver arguments
		obj.odeopt   = struct();
		end % Rietkerk_

		function obj = Rietkerk(varargin)
			obj = obj.Rietkerk_();
			if (~isempty(varargin) && isstruct(varargin{1}))
				obj = copyfields_deep(varargin{1},obj);
			else
			    for idx=1:2:length(varargin)-1
				obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			    end
			end
		end % constructor

	end % methods
end % classdef Rietkerk

