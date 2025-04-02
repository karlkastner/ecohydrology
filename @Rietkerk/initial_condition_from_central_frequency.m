% 2022-01-14 19:34:38.087731800 +0100
%
%% extract dominant frequeny from a previous model run and generate a new
%% initial condition with only this frequency
%% for faster generation of asymptotic patterns
%%
%% for 1D model setups
function [b0,w0,h0] = initial_condition_from_central_frequency(obj,z)
	if (isvector(z))
		z = rvec(z);
	end
	
	[b,w,h]    = obj.extract1(z(end,:));
	
	% determine dominant wave-length
	S          = periodogram_bartlett(b-mean(b),obj.L,round(sqrt(obj.nx)),obj.nx); 
	[Smax,mdx] = max(S);
	fx = fourier_axis(obj.x);
	fc = fx(mdx);

	% round to nearest integer
	fc = round(obj.L*fc)/obj.L;

	if (0) %~obj.opt.legacy_ic)

	[b0,w0,h0] = obj.initial_condition_periodic(fc);

	else

		s = min(mean(b),sqrt(2)*std(b));
		b = mean(b) + s*cos(2*pi*fc*cvec(obj.x));
		% ensure positivity
		b0  = b+min(b);
		w0  = repmat(mean(w),obj.nx,1);
		h0  = repmat(mean(h),obj.nx,1);
	%	z0 = double([b;h;w]);
	end

	% ensure positivity
	%w  = repmat(mean(w),obj.n,1);
	%h  = repmat(mean(h),obj.n,1);
	if (nargout() == 1)
		b0 = double([b0;h0;w0]);
	end
end

