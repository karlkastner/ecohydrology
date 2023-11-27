% Mon  2 May 14:18:38 CEST 2022
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
function [to,zo,out] = solve_split(obj)
	t   = (0:obj.opt.dt:obj.T);
	to  = (0:obj.opt.dto:obj.T);
	nto = length(to);

	% initialize fourier transform of impulse responses of linear
	% advection-diffusion part
	timer = tic();
	obj.init_fourier_matrices();
	out.runtime(1) = toc(timer);
	timer = tic();

	% allocate memory for output
	z  = cvec(obj.opt.compute_class(obj.z0));
	zo = zeros_man(nto,numel(z),obj.opt.output_class);
	zo(1,:) = z; 
	todx = 1;
	for tdx=1:length(t)
		% react half step
		zold = z;
		z = step_react_ralston(t(tdx),0.5*obj.opt.dt,z,@obj.dz_dt_react);
		% advect-diffuse a full step
		z = step_advect_diffuse(t(tdx),z);
		% react half step
		z = step_react_ralston(t(tdx),0.5*obj.opt.dt,z,@obj.dz_dt_react);
		if (t(tdx)>=to(todx+1))
			zo(todx+1,:) = z;
			todx = todx+1;
		end
	end % for tdx
	out.runtime(2) = toc(timer);
	out.y_final = z;

	function z = step_advect_diffuse(t,z)
		zold = z;
		zz = cell(1,obj.nvar);
		[zz{:}] = obj.extract2(z);

		for idx=1:obj.nvar
			if (1 == obj.ndim)
				zz{idx} = ifft(obj.aux.F{idx}.*fft(zz{idx}));
			else
				zz{idx} = flat(ifft2(obj.aux.F{idx}.*fft2(zz{idx})));
			end % if
		end
		z = vertcat(zz{:});
	end	
end % solve_split

