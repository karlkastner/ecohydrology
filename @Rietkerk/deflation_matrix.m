% 2023-07-20 11:37:00.289883536 +0200
% Karl Kastner, Berlin
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
% precompute downsamplin matrices
function [Ax,Ay] = deflation_matrix(obj)
	if (isfield(obj.aux,'Ax_'))
		Ax = obj.aux.Ax_;
		Ay = obj.aux.Ay_;
	else
		Ax = deflation_matrix(obj.n(1))';
		if (obj.ndim>1)
			Ay = downsampling_matrix(obj.n(2),'pairwise');
		else
			Ay = 1;
		end
		obj.aux.Ax = Ax;
		obj.aux.Ay = Ay;
	end
end

