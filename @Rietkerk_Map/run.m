% Mon 31 May 20:20:46 CEST 2021
function [t,y,rk,runtime] = run(obj,varargin)
	% compute hashkey
	rk = Rietkerk(varargin{:});
	key = obj.hash(rk);
	% test if hashkey is not yet in database
	if (~isKey(obj.map,key))
		% add hashkey to database and rewrite table
		obj.map(key) = rk;
		obj.write_table();
	end
	oname = [obj.path_str,filesep,obj.base_str,num2str(key),'.mat'];
	
	if (~exist(oname,'file'))
		printf('Running %d\n',key);
		% run model
		y0 = rk.init();
		timer = tic();
		[t,y] = rk.solve(y0);
		runtime = toc(timer);
		printf('Runtime %g\n',runtime);
		% save disk space
		y = single(y);
		% store model results
		save(oname,'t','y','rk','runtime');
	else
		% load results
		disp(['Loading ',oname])
		load(oname);
	end % else of if ~exist file
end % run
