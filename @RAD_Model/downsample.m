% Tue  7 Nov 10:05:17 CET 2023
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
function obj = resample(obj,direction)
	if (nargin()<2)
		direction = 1;
	end
	% get initial condition
	z0_C = cell(obj.nvar,1);
	[z0_C{:}] = obj.extract2(obj.z0);
	% reset auxiliar variables
	obj.aux = struct();
	nx     = obj.nx;
	f = fieldnames(obj.p);
	% downsampling spatially distributed model coefficients
	if (direction>0)
		if (any(mod(obj.nx,2)))
			error('nx must be even');
		end
		D1 = downsampling_matrix(nx(1));
		D2 = downsampling_matrix(nx(2));
		obj.nx = nx/2;
	else
		% upsampel, the upsampling matrix is just the transpose of the downsampling matrix
		D1 = 2*downsampling_matrix(2*nx(1))';
		D2 = 2*downsampling_matrix(2*nx(2))';
		obj.nx = 2*nx;
	end

	for idx=1:length(f)
		nn = numel(obj.p.(f{idx}));
		if (nn == prod(nx))
			obj.p.(f{idx}) = flat(D1*reshape(obj.p.(f{idx}),nx)*D2');
		end
	end % for idx
	% downsampling the initial condition
	z0_down = [];
	for idx=1:obj.nvar
		z0_down = [z0_down; flat(D1*z0_C{idx}*D2')];
	end % for idx
	obj.z0 = z0_down;
end % downsample

