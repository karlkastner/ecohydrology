% 2021-06-30 13:06:19.787743754 +0200
function b0 = homogeneous_state(obj)
	r = obj.r;
	d = obj.d;
	b0 = [
	                               0,                                              
	 (r + (r^2 - 4*d^2)^(1/2))/(2*d),                                                
	 (r - (r^2 - 4*d^2)^(1/2))/(2*d) ];
end

