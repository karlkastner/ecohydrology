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
function store(obj,t,z)
	% advance to next output slot
	obj.aux.odx = obj.aux.odx+1;
	odx = obj.aux.odx;

	% TODO: zo was allocated, scalars like esum not
	%no = size(obj.out.zo,1);
	no = length(obj.out.esum);
	% reallocate output array
	% only takes effect when output is variable
	% with respect to relative change
	if (odx > no)
		% double vector length this will be truncated at the end of the
		% mode run, if more memory was allocated than needed
		no                  = 2*no;
		% extend vectors by filling zeros
		% to has to be filled with NaN or inf to avoid storing every time step
		obj.out.to(end+1:no) = NaN;
		obj.out.zo(no,1)     = 0;
		obj.out.esum(no,1)     = 0;
		obj.out.n_attempt(no,1)= 0;
		obj.out.n_error_tolerance_exceeded(no,1)= 0;
		obj.out.n_iter(no,1) = 0;
		obj.out.n_liter(no,1) = 0;
		obj.out.n_neg(no,1)  = 0;
		obj.out.n_solver_failed(no,1)= 0;
		obj.out.n_step(no,1)   = 0;
		obj.out.runtime(no,1)  = 0;
	end
	% TODO linear interpolation if it does not match exactly t
	obj.out.to(odx)     = t;
	% since the last flux part can be deactivated, the indices hafe to be specified
	obj.out.zo(odx,1:length(z))   = z;
	% zo_last has to be stored in the same precision as the compute class,
	% and thus has to be stored separately from zo(:,odx)
	obj.aux.zo_last     = z;
	obj.aux.rms_zo_last = rms(obj.aux.zo_last);
	% total runtime
	obj.out.runtime(odx)                    = toc(obj.aux.timer);
end % store

