% 2025-07-28 17:23:55.116572027 +0200
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
function step_postprocess(obj,z)
	tdx = obj.aux.tdx;	
	odx = obj.aux.odx;

	% store time counters summed over time steps over storage intervals
	obj.out.esum(odx) = obj.out.esum(odx) + obj.aux.stat.maxe;
	obj.out.n_attempt(odx) = obj.out.n_attempt(odx) + obj.aux.step.n_attempt;
	obj.out.n_error_tolerance_exceeded(odx) = obj.out.n_error_tolerance_exceeded(odx) + obj.aux.step.n_error_tolerance_exceeded;
	obj.out.n_neg(odx) = obj.out.n_neg(odx) + obj.aux.step.n_neg;
	obj.out.n_iter(odx) = obj.out.n_iter(odx) + sum(obj.aux.step.n_iter);
	if (isfield(obj.aux.step,'n_liter'))
	obj.out.n_liter(odx) = obj.out.n_liter(odx) + sum(obj.aux.step.n_liter);
	end
	obj.out.n_solver_failed(odx) = obj.out.n_solver_failed(odx) + obj.aux.step.n_solver_failed;
	obj.out.n_step(odx) = obj.out.n_step(odx) + 1;
	% 
	if (obj.opt.output.store_step)
		ns = length(out.step);
		if (tdx>ns)
			% reallocate
			obj.out.step(2*ns) = obj.out.step
		end
		obj.out.step(tdx)       = obj.aux.step;
		obj.out.step(tdx).edx   = obj.aux.edx;
		obj.out.step(tdx).t     = t;
		obj.out.step(tdx).dt    = dt;
		obj.out.step(tdx).dt    = obj.aux.dt_max_prec;
		% as this is at the end of step, the values are written at tdx+1
		%obj.out.step(tdx).max = [max(b),max(w),max(h)];
		%obj.out.step(tdx).mean = [mean(b),mean(w),mean(h)];
		%obj.out.step(tdx).min = [min(b),min(w),min(h)];
	end
%	if (obj.opt.store_stat)
%		obj.out.stat(tdx) = stat;
%	end
	if (obj.opt.output.store_state_at_event)
		%[b,w,h]=obj.extract1(z);
		% TODO reallocate
		obj.out.event.z(:,obj.aux.edx) = z;
	end
end

