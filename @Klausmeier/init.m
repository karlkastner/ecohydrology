% 2021-06-30 22:04:16.711272411 +0200
function init(obj)
	n = obj.nx;
	L = obj.L;
	if (size(obj.D1,1)~=n)
		obj.D1 = derivative_matrix_1_1d(n,L,-sign(obj.c),'circular');
		obj.D2 = derivative_matrix_2_1d(n,L,2,'circular');
	end
end
