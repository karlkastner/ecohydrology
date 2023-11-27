% Mon 16 Oct 10:18:57 CEST 2023
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
function [G] = init_fourier_matrices(obj)
	dx  = obj.L./obj.nx;
	% set up advection-diffusion split-solution
	rhs = zeros(prod(obj.nx),1);
	rhs(1,1) = 1;

	% impulse response (fundamental solution)
	G = {};

	% fourier transform of impulse response (transfer function)
	obj.aux.F = {};
	for idx=1:obj.nvar
		% the trapezoidal scheme is only positivity preserving for sufficiently small time steps
		% since the advection diffusion part is linear and time invariant
		% the time step can be split in substeps to satisfy the positivity constraint
		% this does not increase the computationaly effort, as
		% multiple time steps are perfored by taking the element-wise k-th power
		% of the transformation matrix F, which is precomputed
		if (1 == length(obj.nx))
			e = obj.p.ex(idx);
			v = obj.p.vx(idx);
		else
			e = [obj.p.ex(idx),obj.p.ey(idx)];
			v = [obj.p.vx(idx),obj.p.vy(idx)];
		end
		dt_max = 0.5*sum(dx.^2./e);
		k1  = ceil(obj.opt.dt/dt_max);
		k2  = ceil(sqrt(sum((rvec(v).*obj.opt.dt./rvec(dx)).^2)));
		k   = max(k1,k2);
		k   = max(1,k);
		%k(idx) = ceil(obj.opt.dt/dt_max)
		%G{idx} = step_advection_diffusion_trapezoidal(obj.opt.dt/double(k(idx)),rhs,obj.nx,dx,double(v),double(e));
		G{idx} = step_advection_diffusion_trapezoidal(obj.opt.dt/double(k),dx,obj.nx,rhs,double(v),double(e));

		if (1 == length(obj.nx))
			F{idx} = fft(G{idx}).^k;
		else
			F{idx} = fft2(reshape(G{idx},obj.nx)).^k;
		end
		F{idx} = obj.opt.compute_class(F{idx});
		obj.aux.k(idx) = k;
	end
	obj.aux.F = F;
	obj.aux.G = G;
end

