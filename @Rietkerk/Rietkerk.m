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
			  'cb',	10    ... %
			... % range of db given as 0 to 0.5 by rietkerk
			 ,'db',	 0.25 ... %
			 ,'gb',	 0.05 ... %
			 ,'a',	 0.2  ... %
			 ,'ex', {{0.1,0.1,100}} ... % eb,ew,eh
			 ,'ey', {{0.1,0.1,100}} ... % eb,ew,eh
			 ,'vx', {{0,0,0}} ... % vb,vw,vh
			 ,'vy', {{0,0,0}} ... % vb,vw,vh
			 ...% range given as 0 to 3 by Rietkerk
			 ,'R',	 1    ... %
			 ,'w0',	 0.2  ...   % 1
			 ,'kUw',	 5 ...     % mm
			 ,'kIb',	 5 ...     % g m^-2
			 , 'bevap', 0 ...
			 , 'revap', 0.2 ...
			 ...%,'rw',	 0.2  ... % d^-1
			 ...%,'krw', sqrt(eps) ...
			 ...%,'kbrw', 1./sqrt(eps) ...
			 , 'fgb', 1 ... % if not 1, soil water uptake slows down during the dry season
			 , 'kgb', 10 ... % slow down coefficient
			 ... , 'Manning', 0.055 ... % c.f. caviedes 2022 
			 ,'Cb', NaN ...
			 ,'Cv', NaN ...
			 ,'kbC', 0 ...
			 ,'rdrain', 0 ...
			 ,'wdrain', 0 ...
			 ,'lcd', 0 ...
			 ,'zb', 0 ...
			 ,'pI', 1 ...
		); % struct p
		
		obj.unit = struct(  'a', 'd^{-1}' ...
                                  , 'b','g m^{-2}' ...
				  ,'cb','mm^{-1} g m^{-2}' ...
				  ,'db','1/d' ...
				  ,'ex','m^2 d^{-1}' ...
				  ,'gb', 'd^{-1} g^{-1} m^{-2}' ...
				  ,'h','mm' ...
				  ,'R','mm/d' ...
				  ,'vx','m d^{-1}' ...
			  	  ,'w','mm' ...
				  ,'T', 'd' ...
				  ,'L', 'm' ...
				 );

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

		obj.initial_condition = struct( 'mu',[0,0,0] ...
					       ,'sd',[0,0,0] ...
					       ,'sl',[0,0,0]);
		obj.initial_condition.dist = {[],[],[]}; 

		% boundary conditions
		obj.boundary_condition = {'circular','circular'};

		% derivative matrices
		obj.aux      = struct('fgb',1,'fgw',1,'fgh',1);
		obj.aux.surface_flow = true;

		obj.opt.output.base_str = 'rietkerk';
		obj.opt.isreal = true;
		obj.opt.nonlinear_flow = false;
		obj.opt.linear_infiltration = true;
		obj.opt.output.store_fluxes = false;
		obj.opt.output.store_state_at_event = false;
		obj.opt.output.store_stat_at_step = false;
		obj.opt.output.store_step = false;
		obj.opt.event = 'start';

		obj.ptsrel.b = 0;
		obj.ptsrel.w = 0;
		obj.ptsrel.h = 0;
		% TODO move to rad
		obj.opt.temporal_noise = false;

		% solver arguments
		%obj.odeopt   = struct();
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
			if (1 == size(obj.boundary_condition,3))
				for vdx=2:obj.nvar 
					obj.boundary_condition(:,:,vdx) = obj.boundary_condition(:,:,1);
				end % for vdx
			end % if 1 == size(bc)
		end % constructor

	end % methods
end % classdef Rietkerk

