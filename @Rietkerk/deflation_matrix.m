% 2023-07-20 11:37:00.289883536 +0200
% precompute downsamplin matrices
function [Ax,Ay] = deflation_matrix(obj)
	if (isfield(obj.aux,'Ax'))
		Ax = obj.aux.Ax;
		Ay = obj.aux.Ay;
	else
		Ax = deflation_matrix(obj.n(1))';
		if (obj.ndim>1)
			Ay = deflation_matrix(obj.n(2));
		else
			Ay = 1;
		end
		obj.aux.Ax = Ax;
		obj.aux.Ay = Ay;
	end
end

