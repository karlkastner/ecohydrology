% 2021-06-30 22:04:16.711272411 +0200
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
%% this is slow, but only called at startup
function init_advection_diffusion_matrix(obj)
	if (isfield(obj.p,'ex') && ~isempty(obj.p.ex))
	% note, bc of (possible) circular condition, dx is not computee in matrix setup routines
	% as the location of the boundary is interpreted differently
	dx = obj.L./obj.nx;

	% copy boundary from first to other state variables, when only
	% one the bc for the first state variable has been stated

	% stack the matrix of the advection-diffusion part
	nn = prod(obj.nx);

	% identity matrix, first dimension
	obj.aux.Ix   = speye(obj.nx(1));

	% advection diffusion matrix of each state variable, first dimension
	% the splitting schemes operate on each state variable individually
	obj.aux.Axi = {};

	% stacked advection-diffusion matrix, first dimension
	% the alternating implicit direction schmes require dimensionally split matrices
	obj.aux.Ax = spzeros(obj.nvar*nn,obj.nvar*nn);

	if (obj.ndim > 1)
		% identity matrix, second dimension
		obj.aux.Iy   = speye(obj.nx(2));
		% advection diffusion matrix of each state variable, second dimension
		obj.aux.Ayi = {};
		% stacked advection-diffusion matrix, second, dimension
		obj.aux.Ay  = spzeros(obj.nvar*nn,obj.nvar*nn);
	end
	
	% for each state variable
	for vdx=1:obj.nvar
		% TODO, this should go into the matrix setup
if (0)
		% left boundary condition
		bcxl = obj.boundary_condition{1,1,nvar};
		if (iscell(bcxl))
			bcxl = cell2mat(bcxl(1:3));
		end
		% right boundary condition
		bcxr = obj.boundary_condition{1,2,nvar};
		if (iscell(bcxr))
			bcxl = cell2mat(bcxr(1:3));
		end
end
		switch (obj.ndim)
		case {1}
			% first dimension, first derivative (advection)
			% TODO allow for spatially varying advection coefficient
			% TODO combine D1 and D2 setup into one function, to allow for artificial diffusion
			% upwinding
			sign_vx = sign(obj.p.vx{vdx});
			D1x     = -derivative_matrix_1_1d(obj.nx(1),dx(1),sign_vx, ...
					obj.boundary_condition{1,1,vdx}, ...
					obj.boundary_condition{1,2,vdx},true);
			obj.aux.Axi{vdx} = obj.p.vx{vdx}*D1x;
		
			% first dimension, second derivative (diffusion)
			if (isscalar(obj.p.ex{vdx}))
				% constant diffusion
				eD2x = obj.p.ex{vdx}*derivative_matrix_2_1d(obj.nx(1),dx(1),2, ...
					obj.boundary_condition{1,1,vdx},...
					obj.boundary_condition{1,2,vdx},true);
			else
				% spatially varying diffusion coefficient
				eD2x = derivative_matrix_2_1d_vc(obj.nx(1),obj.L(1),obj.p.ex{vdx});
			end
			obj.aux.Axi{vdx} = obj.aux.Axi{vdx} + eD2x;
			% write to stacked matrix
			id = (vdx-1)*nn + (1:nn);
			obj.aux.Ax(id,id) = obj.aux.Axi{vdx};
		case {2} % two dimensions
			%D1y = derivative_matrix_1_1d(obj.nx(2),obj.L(2),-sign(obj.pmu.vx(2)),obj.boundary_condition{1},obj.boundary_condition{2});
			% first dimension, first derivative
			% TODO allow for varying coefficients, pass varying coefficients to the matrix setup directly
			sign_vx = sign(obj.p.vx{vdx});
			vD1x = -obj.p.vx{vdx}*derivative_matrix_1_1d(obj.nx(2),dx(2),sign_vx, ...
							obj.boundary_condition{1,1,vdx}, ...
							obj.boundary_condition{1,2,vdx},true);

			eD2x  = obj.p.ex{vdx}*derivative_matrix_2_1d(obj.nx(1),dx(1),2, ...
							obj.boundary_condition{1,1,vdx}, ...
							obj.boundary_condition{2,2,vdx},true);
			% combined advection diffusion
			obj.aux.Axi{vdx} = kron(vD1x + eD2x, obj.aux.Iy);

			% second dimension, first derivative
			sign_vy = sign(obj.p.vy{vdx});
			vD1y = -obj.p.vy{vdx}*derivative_matrix_1_1d(obj.nx(2),dx(2),sign_vy, ...
							obj.boundary_condition{1,1,vdx}, ...
							obj.boundary_condition{1,2,vdx},true);

			D2y  = obj.o.ey{vdx}*derivative_matrix_2_1d(obj.nx(2),dx(2),2, ...
						obj.boundary_condition{2,1,vdx}, ...
						obj.boundary_condition{2,2,vdx},true);

			% combined advection diffusion
			obj.aux.Ayi{vdx} = kron(vD1y + eD2y, obj.aux.Iy);

			% write to stacked matrix
			id = (vdx-1)*nn + (1:nn);
			obj.aux.Ax(id,id) = obj.aux.Axi{vdx};
			obj.aux.Ay(id,id) = obj.aux.Ayi{vdx};
		end % switch ndim
	end % for vdx

	% combine dimensions
	obj.aux.A = obj.aux.Ax;
	if (obj.ndim > 1)
		obj.aux.A = obj.aux.A + obj.aux.Ay;
	end % if ndim > 1
	end % if ifsfield ex
end % init_advection_diffusion_matrix

