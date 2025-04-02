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
% 2D complex Gizburg-Landau equation on a periodic domain.
% c.f. Kassam & Trefethen
%
% u_t = u - (1+iA)u*abs(u)^2 + D(u_xx + u_yy)
%
classdef Ginzburg_Landau < RAD_Model
	properties
		% all inherited from RAD_Model
	end % properties
	methods
	% pseudo-constructor
	function obj = Ginzburg_Landau_(obj)
		obj.nvar = 1;
		% mean of the physical parameters
		obj.pmu = struct(         ...
			'ex', 1 ...
			,'ey', 1 ...
			,'A', 1.3 ...
			,'vx', 0 ...
			,'vy', 0 ... 
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

		obj.initial_condition = struct( 'mu',[complex(0)] ...
					       ,'sd',[0.1] ...
					       ,'sl',[0] ...	
						);
		obj.initial_condition.dist = {'normal'}; 

		% boundary conditions
		obj.boundary_condition = {'circular','circular'};

		% derivative matrices
		obj.aux      = struct('fgb',1,'fgw',1,'fgh',1);

		obj.opt.base_str = 'ginzburg-landau-';
		obj.opt.isreal = false;

		% solver arguments
		obj.odeopt   = struct();
		end % Ginzburg_Landau_

		function obj = Ginzburg_Landau(varargin)
			obj = obj.Ginzburg_Landau_();
			if (~isempty(varargin) && isstruct(varargin{1}))
				obj = copyfields_deep(varargin{1},obj);
			else
			    for idx=1:2:length(varargin)-1
				obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			    end
			end
		end % constructor

	end % methods
end % classdef Ginzburg_Landau

