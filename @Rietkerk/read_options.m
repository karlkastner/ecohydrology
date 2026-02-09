% 2025-08-03 22:20:03.531618028 +0200
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
function opt = read_options(obj, filename)
	tab = readtable(filename, 'ReadRowNames', true, 'Format', 'auto');
%	opt = rows2vars(tab);
%	opt = table2struct(opt(:,2:end));

	f_C = {
			  'nonlinear_flow', 'double';
			, 'linear_infiltration', 'double';
			, 'loadfinal', 'double';
			, 'adapt_time_step', 'double';
			, 'dt', 'double';
			, 'dt_min', 'double';
			, 'dt_max', 'double';
			, 'outer_abstol', 'double';
			, 'outer_reltol', 'double';
			, 'dt_max_scale_up', 'double';
			, 'dt_min_scale_down', 'double';
			, 'dto', 'double';
			, 'inner_maxiter', 'double';
			, 'compute_class', 'function_handle';
			, 'output_class', 'function_handle';
			, 'path_str', 'char';
			, 'solver', 'char';
			, 'inner_solver', 'char';
			, 'linear_solver', 'char';
			, 'preconditioner', 'char';
			, 'event', 'char';
			, 'store_fluxes', 'double';
			, 'inner_q', 'char';
			, 'zero_inertia.discretization', 'char'
			, 'zero_inertia.input_factor', 'double';
			, 'zero_inertia.output_factor', 'double';
			, 'rng', 'double';
			, 'dt_event', 'double';
			, 'time_stepper', 'char';
	};
	opt = struct()

	for idx=1:height(tab)
		rowname = tab.Row{idx};
		if (rowname(1) ~= '#')
		val = tab{idx,1}{1};
		id = find(arrayfun(@(x) strcmp(rowname,x), f_C(:,1)));
		switch (length(id))
		case {0}
			error(['option "', rowname, '" is not defined'])
		case {1}
			switch (f_C{id,2})
			case ('double')
				val = str2num(val);
			case ('function_handle')
				val = str2func(val);
			case ('char')
				% nothing to convert
			end
			opt = setfield_deep(opt,rowname,val);	
		otherwise
			error(['option', rowname, 'is defined more than once'])
		end
		end
	end
end

