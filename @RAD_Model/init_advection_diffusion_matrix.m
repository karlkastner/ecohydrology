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
function init_advection_diffusion_matrix(obj)
	if (isfield(obj.p,'ex') && ~isempty(obj.p.ex))
	% note, bc of circular condition, do not compute in matrix setup routines
	dx = obj.L./obj.nx;

	bcl = obj.boundary_condition{1,1};
	if (iscell(bcl))
		bcl = cell2mat(bcl(1:3));
	end

	% first dimension
	obj.aux.D1xl = derivative_matrix_1_1d(obj.nx(1),dx(1),+1,bcl,obj.boundary_condition{1,2},true);
	obj.aux.D1xr = derivative_matrix_1_1d(obj.nx(1),dx(1),-1,bcl,obj.boundary_condition{1,2},true);
	obj.aux.D2x  = derivative_matrix_2_1d(obj.nx(1),dx(1),2,bcl,obj.boundary_condition{1,2},true);
	obj.aux.D1xl1 = obj.aux.D1xl;
	obj.aux.D1xr1 = obj.aux.D1xr;
	obj.aux.D2x1  = obj.aux.D2x;

	% second dimension
	if (obj.ndim > 1)
		%D1y = derivative_matrix_1_1d(obj.nx(2),obj.L(2),-sign(obj.pmu.vx(2)),obj.boundary_condition{1},obj.boundary_condition{2});
		obj.aux.D1yl1 = derivative_matrix_1_1d(obj.nx(2),dx(2),+1,obj.boundary_condition{2,1},obj.boundary_condition{2,2},true);
		obj.aux.D1yr1 = derivative_matrix_1_1d(obj.nx(2),dx(2),-1,obj.boundary_condition{2,1},obj.boundary_condition{2,2},true);
		obj.aux.D2y1  = derivative_matrix_2_1d(obj.nx(2),dx(2),2,obj.boundary_condition{2,1},obj.boundary_condition{2,2},true);
		obj.aux.Ix   = speye(obj.nx(1));
		obj.aux.Iy   = speye(obj.nx(2));
		obj.aux.D1xr = kron(obj.aux.D1xr1,obj.aux.Iy);
		obj.aux.D1xl = kron(obj.aux.D1xl1,obj.aux.Iy);
		obj.aux.D1yr = kron(obj.aux.Ix,obj.aux.D1yr1);
		obj.aux.D1yl = kron(obj.aux.Ix,obj.aux.D1yl1);
		obj.aux.D2x  = kron(obj.aux.D2x1,obj.aux.Iy);
		obj.aux.D2y  = kron(obj.aux.Ix,obj.aux.D2y1);
	end
	% stack the matrix of the advection-diffusion part
	n = prod(obj.nx);
	%obj.aux.AD = spzeros(obj.nvar*n,obj.nvar*n);
	obj.aux.AD = [];
	obj.aux.Ax = [];
	obj.aux.Ay = [];
	% spalloc(obj.nvar*n,obj.nvar*n,5*obj.nvar*n);
	Z = spzeros(n,n);
	% matrix entries for each state variable
	for idx=1:obj.nvar
		% first dimension
		% diffusion part
		ADxi  = obj.p.ex(idx)*obj.aux.D2x;
		ADx1i = obj.p.ex(idx)*obj.aux.D2x1;
		% advection part
		if (obj.p.vx(idx) < 0)
			ADxi  = ADxi  + obj.p.vx(idx)*obj.aux.D1xl;
			ADx1i = ADx1i + obj.p.vx(idx)*obj.aux.D1xl1;
		else
			ADxi  = ADxi  + obj.p.vx(idx)*obj.aux.D1xr;
			ADx1i = ADx1i + obj.p.vx(idx)*obj.aux.D1xr1;
		end
		if (length(obj.nx)>1)
			% diffusion part
			ADyi  = obj.p.ey(idx)*obj.aux.D2y;
			ADy1i = obj.p.ey(idx)*obj.aux.D2y1;
			% advection part		
			if (obj.p.vy(idx) < 0)
				ADyi  = ADyi  + obj.p.vy(idx)*obj.aux.D1yl;
				ADy1i = ADy1i + obj.p.vy(idx)*obj.aux.D1yl1;
			else
				ADyi  = ADyi  + obj.p.vy(idx)*obj.aux.D1yr;
				ADy1i = ADy1i + obj.p.vy(idx)*obj.aux.D1yr1;
			end	
		end
		% write to matrix comprising of all dimensions
		% this is slow, but only called at startup
		% TODO this can be more elegantly written with kron([0,1,0],A)
		% obj.aux.AD((idx-1)*n+1:idx*n,(idx-1)*n+1:idx*n) = ADi; 
		A = [];
		Ax = [];
		Ay = [];
		for jdx=1:idx-1
			A  = [A,Z];
			Ax = [Ax,Z];
			if (obj.ndim>1)
				Ay = [Ay,Z];
			end
		end
		Ax = [A,ADxi];
		if (obj.ndim > 1)
			A = [A,ADxi+ADyi];
			Ay = [Ay,ADyi];
		else
			A = [A,ADxi];
		end
		for jdx=idx+1:obj.nvar
			A = [A,Z];
			Ax = [Ax,Z];
			if (obj.ndim>1)
				Ay = [Ay,Z];
			end
		end
		obj.aux.AD = [obj.aux.AD;A];
		obj.aux.Ax = [obj.aux.Ax;Ax];
		%obj.aux.Ax{idx} = ADxi;
		obj.aux.Ax1{idx} = ADx1i;
		if (obj.ndim > 1)
			obj.aux.Ay = [obj.aux.Ay;Ay];
			% obj.aux.Ay{idx} = ADyi;
			obj.aux.Ay1{idx} = ADy1i;
		end
	end
	end % if ~isepmpty ex
	obj.aux.I1 = speye(prod(obj.nx));
	obj.aux.I = speye(obj.nvar*prod(obj.nx));
end % init_advection_diffusion_matrix

