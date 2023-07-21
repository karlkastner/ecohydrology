% Thu 20 Jul 10:20:20 CEST 2023
% downsample all 2D variables with factor of 1/2
function obj2 = downsample(obj)
	obj2   = Rietkerk(obj);
	if (max(mod(obj.n,2))>0)
		error('can only downsample with even number of grid points');
	end
	obj2.n = obj.n/2;
	% initial conidition
	obj2.z0 = obj.downsample_z(obj.z0);
	field_C=fieldnames(obj.p);
	[Ax,Ay] = obj.deflation_matrix();
	for idx=1:length(field_C);
		siz = size(obj.p.(field_C{idx}));
		if (siz(1) == obj.n(1) && (obj.ndim == 1 || siz(2) == obj.n(2)))
			obj2.p.(field_C{idx}) = Ax*(obj.p.(field_C{idx}))*Ay;
		end
	end
	% reset aux
	obj2.aux = struct();
end


