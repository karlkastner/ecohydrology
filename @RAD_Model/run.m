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
function [t,y,out] = run(obj) %varargin)
	% unqiue identifier for parameter combinations
	key = obj.hash();

	[oname,oname_final] = obj.filename();
	if (  ~exist(oname_final,'file') )
		printf('Running %d\n',key);
		printf([oname_final,'\n']);
		% reserve output files (quasi-semapthore for parallel computation)
		y = [];
		t = [];
		runtime = [];
		out = struct('y', y,'t',t,'y_final',[],'runtime',[]);
		obj.save(t,y,out);
		% run model
		tic();
		obj.init();
		runtime_init = toc();
		[t,y,out]   = obj.solve();
		out.runtime = [runtime_init,out.runtime];
		printf('Runtime %g\n',out.runtime(end));
		obj.p
		obj.save(t,y,out);
	else
		% load results
		printf('Loading %d\n',key);
		[t,y,out] = obj.load();
		printf('Runtime %g\n',out.runtime(end));
	end % else of if ~exist file
end % run

