% 2022-01-14 19:34:38.087731800 +0100
%
%% extract dominant frequeny from a previous model run and generate a new
%% initial condition with only this frequency
%% for faster generation of asymptotic patterns
function y0 = initial_condition_from_central_frequency(obj,y)
	if (isvector(y))
		y = rvec(y);
	end
	[b,w,h]    = obj.extract1(y(end,:));
	% determine dominant wave-length
	S          = periodogram_bartlett(b-mean(b),obj.L,round(sqrt(obj.n)),obj.n); 
	[Smax,mdx] = max(S);
	fx = fourier_axis(obj.x);
	fc = fx(mdx);
	% round to nearest integer
	fc = round(obj.L*fc)/obj.L;
	%b  = 0.01*mean(b).*0.5*(1 + cos(2*pi*fc*cvec(obj.x)));
	% limit variance to ensure positivity of b
	s = min(mean(b),sqrt(2)*std(b));
	b = mean(b) + s*cos(2*pi*fc*cvec(obj.x));
	% ensure positivity
	b  = b+min(b);
	w  = repmat(mean(w),obj.n,1);
	h  = repmat(mean(h),obj.n,1);
	y0 = double([b;h;w]);
end

