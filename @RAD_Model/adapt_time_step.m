% Thu  7 Dec 16:41:17 CET 2023
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
function dt = adapt_time_step(obj,t,dt_old,dt0,stat)
	% max time step satisfying the error constraint
	dt = dt0;

 	% limit decrease
	dt = max(dt,dt_old*obj.opt.dt_min_scale_down);

	% limit increase
	dt = min(dt,dt_old*obj.opt.dt_max_scale_up);

	% limit at lower bound
	dt = max(dt,obj.opt.dt_min);

	% limit at upper bound
	dt = min(dt,obj.opt.dt_max);

	% avoid too many inner steps
	switch (obj.opt.inner_solver)
	case {'newton'}
	switch (obj.opt.innersolver2)
	case {'multigrid-java'}
		iter_max = obj.opt.inner2_maxiter_;
		if (~isempty(iter_max))
			% the last one is the hardest always
			iter = stat.inner2.iter(end);
			dt = min(dt,dt_old*iter_max/iter);
		end
	end % switch i2
	end % switch innersolver

	% limit final step to final time
	dt = min(dt,obj.T-t);
end

