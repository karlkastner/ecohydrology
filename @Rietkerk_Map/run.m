% Mon 31 May 20:20:46 CEST 2021
% Karl Kästner, Berlin
%
%% run the Rietkerk model with parameters specified by varargin,
%% or retrieve the saved results, when the model was already run
function [t,y,rk,runtime] = run(obj,varargin)
	% compute hashkey
	rk  = Rietkerk(varargin{:});
	key = obj.hash(rk);
	if (0)
	% test if hashkey is not yet in databas
	if (~isKey(obj.map,key))
		% add hashkey to database and rewrite table
		obj.map(key) = rk;
		obj.write_table();
	end
	end
	oname       = [obj.path_str,filesep,obj.base_str,num2str(key),'.mat'];
	oname_final = [obj.path_str,filesep,obj.base_str,num2str(key),'-final.mat'];
	
	if (   ~exist(oname,'file') ...
            && ~(obj.loadfinal && exist(oname_final,'file')))
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
		save(oname,'-v7.3','t','y','rk','runtime');
		y_ = y;
		t_ = t;
		y = y([1,end],:);
		t = t([1,end]);
		save(oname_final,'-v7.3','t','y','rk','runtime');
		y = y_;
		t = t_;
	else
		% load results
		if (~obj.loadfinal)
			disp(['Loading ',oname]);
			load(oname);
		else
			disp(['Loading ',oname_final]);
			load(oname_final);
			% TODO hot fix
			if (size(t,1) == 2 && size(t,2)>1)
				t = t(1,[1,end])';
			end
		end
		printf('Runtime %g\n',runtime);
	end % else of if ~exist file
end % run
