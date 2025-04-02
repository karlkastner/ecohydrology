% Wed 18 Oct 08:42:39 CEST 2023
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
classdef Dodorico < RAD_Model
	properties
	end
	methods
		function obj = Dodorico(varargin)
			obj = obj.Dodorico_();
			%RAD_Model@RAD_Model(obj)
			if (~isempty(varargin) && isstruct(varargin{1}))
				obj = copyfields_deep(varargin{1},obj);
			else
			    for idx=1:2:length(varargin)-1
				obj = setfield_deep(obj,varargin{idx},varargin{idx+1});
			    end
			end
		end
		function obj = Dodorico_(obj)
			obj.nvar = 1;
			obj.pmu = struct();
			obj.pmu.ze   = 0.0001;
			obj.pmu.zmax = 1; 
			obj.pmu.ex   = 0.3;
			obj.pmu.ey   = 0.3;
			obj.pmu.vx   = 0.0;
			obj.pmu.vy   = 0.0;
			obj.pmu.a    = 0.45;
			obj.pmu.f    = 0.00;
			obj.pmu.g    = 0;
			obj.pmu.l0   = 0.65;	
			obj.pmu.b    = -0.9*0.65;
			obj.pmu.w0   = 0.04;

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
			end % for 
			% obj.psdist.w0   = 'exp';

			obj.initial_condition = 'random';

			obj.opt           = struct();
			obj.opt.rng       = 0;
			obj.opt.solver    = @ode23;
			obj.opt.path_str  = './';
			obj.opt.base_str  = 'dodorico-';
			obj.opt.loadfinal = false;
			opt.mode          = 'continuous';

			obj.boundary_condition       = {'circular','circular'};

		end % Dodorico_
		function obj = init(obj)
			obj.init@RAD_Model();
			switch (obj.opt.mode)
			case {'discrete'}
				obj.aux.next = zeros(prod(obj.nx),1);
			end
		end
	end % methods
	methods (Static)
		[z0,J] = derive_homogeneous_state();
		[flag,z0,dz_dt] = test_homogeneous_state();
		zp = derive_poles();
	end
end

