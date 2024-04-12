% 2024-02-20 18:00:14.707860259 +0100
function id = reorder(obj,nn)
	nn = prod(obj.nx)
	id=(1:nn)';
	id=id+(0:obj.nvar-1)*nn;
	id=flat(id');
end

