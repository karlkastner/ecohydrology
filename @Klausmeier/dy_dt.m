% 2021-06-30 22:04:16.711272411 +0200
function dy_dt = klausmeier_dt(obj,t,x,y,c,r,d)
	[b,w] = obj.extract(y);

	if (0)
		r = r*(1 + e*randn(obj.n,size(b,2)));
	end
	if (0)
		%d = d*(1 + e*randn(obj.n,size(b,2)));
		d = d*gamrnd(1/e,e,obj.n,size(b,2));
	end
	if (0)
		%u = u*gamrnd(1/e,e,n,1);
		u = u*gamrnd(1/e^2,e^2,n,1);
	end
	if (0)
		v = v*(1 + e*randn(obj.n,size(b,2)));
	end
	if (e>0)
		eb = e*randn(n,1).*b;
	else
		eb = 0;
	end

	db_dt = obj.u.*w.*b.^2 - obj.d.*b   + obj.D2*b + eb;
	dw_dt = obj.r - w - obj.u.*w.*b.^2 + obj.v.*(obj.D1*w);

	if (0)
		db_dt = db_dt+e*randn(obj.n,size(b,2));
	end

	dy_dt = [db_dt; dw_dt];
end
