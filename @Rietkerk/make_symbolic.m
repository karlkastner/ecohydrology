% 2021-07-02 13:29:01.973628215 +0200
% Karl Kästner, Berlin
%
%% make model parameters symbolic
%
function [p,s] = make_symbolic(obj)
	f_C = fieldnames(obj.pmu);
	p = struct();
	for idx=1:length(f_C)
		p.(f_C{idx}) = sym(f_C{idx});
	end
	obj.p = p;
	obj.D1c = sym('D1x');
	obj.D1c = sym('D1x');
	obj.D2  =  sym('D2x');
	obj.D1  = sym('D1');
	obj.I   = 1;
	obj.Z   = 0;
end
