% Mon 16 Oct 09:40:37 CEST 2023
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
function [t,y,out] = solve_euler_forward()
	out = struct('runtime',[]);
	tic();
	obj.init_advection_diffusion_matrix();
	out.runtime(1) = toc();
	tic();

	t = (0:obj.opt.dt:obj.T);
	nt = length(t);
	zz = zeros(prod(obj.nx),nt,obj.opt.output_class);
	zz(:,1) = obj.z0;
	for idx=2:nt
		dz_dt = obj.dz_dt(t(idx-1),z);
		z = z + obj.opt.dt*dz_dt;
		zz(:,idx) = z;
	end
	out.y_final = z;
	zz = zz.';
	out.runtime(2) = runtime;
end

