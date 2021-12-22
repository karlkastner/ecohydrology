% 2021-07-02 13:29:01.973628215 +0200

% make parameters symbolic
function [p,s] = make_symbolic(obj)
	%'cb','db','kb','gb','eb','a','rw','ew','Dh','eh','R','w0','k'};
	f_C = fieldnames(obj.pmu);
	p = struct();
	for idx=1:length(f_C)
		p.(f_C{idx}) = sym(f_C{idx});
	end
	obj.p = p;
	obj.D1c = sym('D1x');
	obj.D1c = sym('D1x');
	obj.D2  =  sym('D2x');
end
