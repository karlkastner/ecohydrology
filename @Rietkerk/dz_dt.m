% Mon 31 May 20:20:46 CEST 2021
%
% p : parameter vector
% s : standard deviation of paramter
% TODO variation of dH is not well implemented good -> define at interfaces between cells
function dz_dt = dz_dt(obj,t,z)
	p = obj.p;
	s = obj.pst;
	if (size(z,2)>1)
	[b,w,h] = obj.extract2(z);
	else
	[b,w,h] = obj.extract1(z);
	end
	if (~isvector(z))
		b=b';
		w=w';
		h=h';
	end
		
	n = prod(obj.n);

		dt = 1; dx=1;
		if (s.db > 0)
			db = p.db*gamrnd(1/(s.db*dt*dx),s.db*(dt*dx),length(x),1);
		else
			db = p.db;
		end
		if (s.eb > 0)
			eb  = p.eb*gamrnd(1/(s.eb*dt*dx),s.eb*dx*dt,length(x),1);
		else
			eb = p.eb*ones(n,1);
		end
		if (s.vh > 0)
			vh  = p.vh*gamrnd(1/(s.vh*dt*dx),s.vh*dx*dt,length(x),1);
		else
			vh = p.vh.*ones(n,1);
		end
		if (s.cb > 0)
			cb  = p.cb*gamrnd(1/(s.cb*dt*dx),s.cb*dx*dt,length(x),1);
		else
			cb = p.cb;
		end
		if (s.a > 0)
			a  = p.a.*gamrnd(1/(s.a*dt*dx),(s.a*dx)*dt,length(x),1);
		else
			a = p.a;
		end
		if (s.rw > 0)
			rw  = p.rw.*gamrnd(1/(s.rw*dt*dx),(s.rw*dx)*dt,length(x),1);
		else
			rw = p.rw;
		end
		if (s.gb > 0)
			gb  = p.gb*gamrnd(1/(s.gb*dt*dx),(s.gb*dx)*dt,length(x),1);
		else
			gb = p.gb;
		end
		if (s.R > 0)
			R  = p.R*gamrnd(1/(s.R*dt*dx),(s.R*dx)*dt,length(x),1);
		else
			R = p.R;
		end
	
	% uptake of water by plants
	U = gb.*w./(w + p.kw).*b;

	% infiltration of water into soil
%	I = a.*h.*(b + p.kb.*p.w0)./(b+p.kb);
	I = a.*obj.infiltration_enhancement(b).*h;

	db_dt = cb.*U - db.*b + eb.*(obj.D2*b) + (obj.D1c*eb).*(obj.D1c*b);
	%db_dt = cb*U - db_.*b + eb_.*(D2*b);
	dw_dt = I - U - rw.*w + p.ew.*(obj.D2*w);
	dh_dt = p.R - I + vh.*(obj.D1*h) + p.eh*obj.D2*h;
	% stack output
	dz_dt = [db_dt; dw_dt; dh_dt];
end % dz_dt

