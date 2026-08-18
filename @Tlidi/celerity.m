% Thu 11 Jun 19:03:05 CEST 2026
function c = celerity(obj,t,b)
	dx = obj.dx;
	db_dt = obj.dz_dt(t,b(:));
	if (obj.ndim>1)
		db_dt = reshape(db_dt,obj.nx);
	end
	db_dx = cdiff(b)/dx(1);
	
	%s = quantile(abs(db_dx),0.95);
	%fdx = (db_dx>s);
	% dz_dt = a*dz_dx
	% a = dz_dt/dz_dx
	%fdx = flat(b>0.1 & b<0.8);
	%c = db_dx(fdx) \ dz_dt(fdx);
	c = db_dx(:) \ db_dt(:);
end

