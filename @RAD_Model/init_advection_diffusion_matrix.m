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
	% note, bc of circular condition, do not compute in matrix setup routines
	dx = obj.L./obj.nx;

	% first dimension
	obj.aux.D1xl = derivative_matrix_1_1d(obj.nx(1),dx(1),+1,obj.boundary_condition{1},obj.boundary_condition{2},true);
	obj.aux.D1xr = derivative_matrix_1_1d(obj.nx(1),dx(1),-1,obj.boundary_condition{1},obj.boundary_condition{2},true);
	obj.aux.D2x  = derivative_matrix_2_1d(obj.nx(1),dx(1),2,obj.boundary_condition{1},obj.boundary_condition{2},true);

	% second dimension
	if (length(obj.nx) > 1)
		%D1y = derivative_matrix_1_1d(obj.nx(2),obj.L(2),-sign(obj.pmu.vx(2)),obj.boundary_condition{1},obj.boundary_condition{2});
		obj.aux.D1yl = derivative_matrix_1_1d(obj.nx(2),dx(2),+1,obj.boundary_condition{1},obj.boundary_condition{2},true);
		obj.aux.D1yr = derivative_matrix_1_1d(obj.nx(2),dx(2),-1,obj.boundary_condition{1},obj.boundary_condition{2},true);
		obj.aux.D2y = derivative_matrix_2_1d(obj.nx(2),dx(2),2,obj.boundary_condition{1},obj.boundary_condition{2},true);
		Ix = speye(obj.nx(1));
		Iy = speye(obj.nx(2));
		obj.aux.D1xr = kron(obj.aux.D1xr,Iy);
		obj.aux.D1xl = kron(obj.aux.D1xl,Iy);
		obj.aux.D1yr = kron(Ix,obj.aux.D1yr);
		obj.aux.D1yl = kron(Ix,obj.aux.D1yl);
		obj.aux.D2x = kron(obj.aux.D2x,Iy);
		obj.aux.D2y = kron(Ix,obj.aux.D2y);
	end
	% stack the matrix of the advection-diffusion part
	n = prod(obj.nx);
	%obj.aux.AD = spzeros(obj.nvar*n,obj.nvar*n);
	obj.aux.AD = [];
	% spalloc(obj.nvar*n,obj.nvar*n,5*obj.nvar*n);
	Z = spzeros(n,n);
	% matrix entries for each state variable
	for idx=1:obj.nvar
		% first dimension
		% diffusion part
		ADi = obj.p.ex(idx)*obj.aux.D2x;
		% advection part
		if (obj.p.vx(idx) < 0)
			ADi = ADi + obj.p.vx(idx)*obj.aux.D1xl;
		else
			ADi = ADi + obj.p.vx(idx)*obj.aux.D1xr;
		end
		if (length(obj.nx)>1)
			% diffusion part
			ADi = ADi + obj.p.ey(idx)*obj.aux.D2y;
			% advection part		
			if (obj.p.vy(idx) < 0)
				ADi = ADi + obj.p.vy(idx)*obj.aux.D1yl;
			else
				ADi = ADi + obj.p.vy(idx)*obj.aux.D1yr;
			end	
		end
		% write to matrix comprising of all dimensions
		% this is slow:
		% obj.aux.AD((idx-1)*n+1:idx*n,(idx-1)*n+1:idx*n) = ADi; 
		A = [];
		for jdx=1:idx-1
			A = [A,Z];
		end
		A = [A,ADi];
		for jdx=idx+1:obj.nvar
			A = [A,Z];
		end
		obj.aux.AD = [obj.aux.AD;A];
	end
	obj.aux.I = speye(obj.nvar*prod(obj.nx));
end

