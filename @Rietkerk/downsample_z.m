% 2023-07-20 11:38:04.779000680 +0200
% downsample single variable by factor 1/2
% function z2 = downsample_z(obj,z)
function z2 = downsample_z(obj,z)
	[Ax,Ay] = obj.deflation_matrix();
	[b,w,h] = obj.extract2(z);
	if (isvector(b))
		b = (Ax*b);
		w = (Ax*w);
		h = (Ax*h);
	else
		b = Ax*(b*Ay);
		w = Ax*(w*Ay);
		h = Ax*(h*Ay);
	end
	z2 = [flat(b);flat(w); flat(h)];
end


