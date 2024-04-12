% Tue 17 Oct 09:39:29 CEST 2023
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
% template class for Reaction-Advection-Diffusion models
classdef RAD_Model < handle
	properties
		% number of grid points
		nx
		% number of variables
		nvar
		% side length of computational domain
		L
		% duration of model run
		T
		% number of output time steps
		%dt
		% output time step
		%dto
		% physical parameters of the Rietkerk model
		p   = struct();
		% mean of physical paramters
		pmu
		% spatial standard deviation of parameters
		pss
		% spatial correlation length of parameters
		psl
		% distribution type
		psdist

		initial_condition
		% boundary condition
		boundary_condition = {'circular','circular'};
		% initial condition
		z0	
		% options
		opt = struct('heterogeneity_integration_order',9 ...
		             ,'rng', 0 ...
		             ,'dt_min', 0 ...
		             ,'dt_max', inf ...
			     ,'dt_scale', sqrt(2) ...
			     ,'inner_m', 10 ...
			     , 'inner_tol', 1e-3 ...
			     , 'inner_solver', 'gauss-newton' ...
			     , 'innersolver2', 'bicgstabl' ...
			     , 'inner_q', 0.5 ...
			     , 'inner_maxiter', 10 ...
			     , 'inner2_maxiter', 1e3 ...
			     , 'inner2_tol', sqrt(eps) ...
			     , 'inner2_maxiter_', inf ...
			     , 'compute_class', @single ...
			     , 'output_class', @half ...
			     , 'loadfinal', true ...
			     , 'outer_abstol', 1e-3 ...
			     , 'outer_reltol', 1e-3 ...
			     , 'solver', 'solve_split' ...
			     , 'preconditioner', 'multigrid-java' ...
			     , 'extrapolation_order', 0 ...
			);
		% solver options
		odeopt
		% fields to be used to create unique model hash
		hashfield_C
		aux
		out
	end
	methods
		function obj = RAD_Model(varargin)
			if (~isempty(varargin) && isstruct(varargin{1}))
				obj = copyfields_deep(varargin{1},obj);
			else
			    for idx=1:2:length(varargin)-1
				obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			    end
			end
		end

		function [x,y] = x(obj,id)
			if (nargin()<2)
				id = 1;
			end
			x = (0:obj.nx(id)-1)*obj.L(id)/obj.nx(id);
			if (nargout()>1)
				y =  (0:obj.nx(2)-1)*obj.L(2)/obj.nx(2);
			end
			%if (length(obj.nx)>1)
			%	y = (0:obj.nx(2)-1)*obj.L(2)/obj.nx(2);
			%end
		end
		function dx = dx(obj)
			dx = rvec(obj.L)./rvec(obj.nx);
		end
		function [fx,fy,frr] = fourier_axis(obj)
			if (obj.ndim == 1)
				fx = fourier_axis(obj.L,obj.nx);
			else
				[fx,fy,frr] = fourier_axis_2d(obj.L,obj.nx);
			end
		end
		function ndim = ndim(obj)
			ndim = length(obj.nx);
		end
	end
end

