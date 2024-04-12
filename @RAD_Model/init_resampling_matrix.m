function init_resampling_matrix(obj)
		[obj.aux.R1down, obj.aux.R1up] = downsampling_matrix_2d(obj.nx);
		if (any(mod(obj.nx,2)))
			error('nx must be even');
		end

end

