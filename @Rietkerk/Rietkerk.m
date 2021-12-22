% Mon  5 Jul 18:00:15 CEST 2021
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
% c.f. Rietkerk et al. 2002, Self-Organization of Vegetation in Arid Ecosystems 
%
classdef Rietkerk < handle
	properties
		% number of grid points
		n
		% transect length
		L
		% duration of model run
		T
		% number of output time steps
		nt
		% output time step
		%dt
		% physical parameters of the Rietkerk model
		p   = struct();
		pmu = struct(         ...
			  'cb',	10    ... % g mm^-1 m^-2
			... % range given as 0 to 0.5 by rietkerk
			 ,'db',	 0.25 ... % d^-1
			 ,'kw',	 5    ... % mm
			 ,'gb',	 0.05 ... % mm g^-1 m^-2 d^-1
			 ,'eb',	 0.1  ... % m^2 d^-1
			 ,'a',	 0.2  ... % d^-1
			 ,'rw',	 0.2  ... % d^-1
			 ,'ew',	 0.1  ... % m^2 d^-1
			 ,'vh',	 0.0  ... %  
		         ... % eh is 100 in the rietkerk model
		 	 ,'eh',	 100    ... % m^2 d^-1
			 ...% range given as 0 to 3 by Rietkerk
			 ...% R was initiall 1.1
			 ,'R',	 1    ... % mm d^-1
			 ,'w0',	 0.2  ...   % 1
			 ... % TODO used to be 1 in the first runs
			 ,'kb',	 5 ...     % g m^-2
		); % struct p
		% standard deviation of parameters per unit time
		% (continuous perturbation)
		pst    = struct(   ...
			  'cb',	0 ...
			 ,'db',	0 ...
			 ,'kw', 0 ...
			 ,'gb',	0 ...
			 ,'eb',	0 ... 
			 ,'a',	0 ... 
			 ,'rw',	0 ... 
			 ,'ew',	0 ... 
			 ,'vh',	0 ...
			 ,'eh',	0 ...
			 ,'w0',	0 ...
			 ,'R',  0 ... 
			 ,'kb', 0 ... 
		); % pst
		% standard deviation of parameters per unit distance
		% (perturbation once during initialization)
		pss    = struct(   ...
			  'cb',	0 ...
			 ,'db',	0 ...
			 ,'kw', 0 ...
			 ,'gb',	0 ...
			 ,'eb',	0 ... 
			 ,'a',	0 ... 
			 ,'rw',	0 ... 
			 ,'ew',	0 ... 
			 ,'vh',	0 ...
			 ,'eh',	0 ...
			 ,'w0',	0 ...
			 ,'R',  0 ... 
			 ,'kb', 0 ... 
		); % pss
		initial_condition = 'random';
		% state of random number generator at start
		rng = 0;
		% derivative matrices
		D1
		D1c
		D2
		% identity matrix
		I
		% zero matrix
		Z
		% time stepping algorithm
		solver   = @ode23;
		solver_stationary = @picard;
		odeopt   = struct();
		symbolic = false;
		order    = 2;
	end % properties
	methods
		function obj = Rietkerk(varargin)
			for idx=1:2:length(varargin)-1
				obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			end
		end % constructor
		function x = x(obj,id)
			if (nargin()<2)
				id = 1;
			end
			x = (0:obj.n(id)-1)*obj.L(id)/obj.n(id);
			%if (length(obj.n)>1)
			%	y = (0:obj.n(2)-1)*obj.L(2)/obj.n(2);
			%end
		end
		function dx = dx(obj)
			dx = obj.L/obj.n;
		end
		function dt = dt(obj)
			dt = obj.T/obj.nt;
		end
	end % methods
end % classdef Rietkerk

