% Mon 31 May 20:20:46 CEST 2021
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
%
%% run the Rietkerk model with parameters specified by varargin,
%% or retrieve the saved results, when the model was already run
function [t,z,out] = run(obj)
	obj.check_parameters();

	% unqiue identifier for parameter combinations
	%key = obj.hash();
	[oname,oname_final,key] = obj.filename();
	if (  ~exist(oname_final,'file') )
%		try
			printf('Running %d\n',key);
			printf([oname_final,'\n']);
			% reserve output files (quasi-semapthore for parallel computation)
			% store aux
			aux = obj.aux;
			obj.save();
			obj.aux = aux;
			% run model
			tic();
			obj.init();
			obj.out.runtime_init = toc();
			[t,z,out]   = obj.solve();
			printf('Runtime %g\n',obj.out.runtime_init+obj.out.runtime(end));
			out = obj.save();
			if (obj.opt.loadfinal)
				z = z([1,end],:)';
			end
%		catch exception
%			% store the exception
%			save([oname(1:end-4),'-exception.mat'],'exception');
%			rethrow(exception)
%		end
	else
		% load results
		if (obj.opt.load_existing)
			printf('Loading %d\n',key);
			[t,z,out] = obj.load();
			if (isfield(out,'runtime') && ~isempty(out.runtime))
				printf('Runtime %g\n',out.runtime(end));
			end
		else
			t = [];
			z = [];
			out = struct();
			printf('Output file %d exists, not running\n',key);
		end
	end % else of if ~exist file
end % run

