% Tue  3 May 13:12:57 CEST 2022
if (0)		
			param = struct();
			param.L = 1e3;
			param.n = 1e3;
			param.T =1000;
			param.nt = 2;



			rk = Rietkerk(param);
			y0 = rk.init();
			tic
%			[t,y] = rk.solve(y0);
%			y_ref = y;
			length(t)
			toc

			figure(1)
			clf
			subplot(2,2,1)
			plot(y([end],:)')
			hold on


			dt = [2,1,0.5,0.25,0.125,0.0625];
			dt = fliplr(dt);
			r = [];
			for idx=1:length(dt)
				t = linspace(0,param.T,round(param.T)/dt(idx));
				tic()
				[t,y] = rk.solve_split(y0,t);
				rt(idx) = toc;
				plot(y(end,:))
				if (1==idx)
					y1 = y;
				end
				r(idx,1) = rms(y(end,:) - y_ref(end,:))	
				r(idx,2) = rms(y(end,:) - y1(end,:))	
			end
			r = r/rms(y_ref(end,:));
			subplot(2,2,2)
			loglog(dt,r,'.-')
			% r = c(1)*dt^c(2)
			c=[ones(length(dt)-1,1),log(cvec(dt(2:end)))] \ log(r(2:end,2)); hold on;
			plot(dt,exp(c(1)).*dt.^c(2))
			subplot(2,2,3)
			lolog(dt,rt)
end

			%t = linspace(t(1),t(end),2*length(t));
			%[t,y] = rk.solve_split(y0,t);
			%toc
			%plot(y(end,:))	
%		subplot(2,2,2)
%		plot(t,1:length(t))
		

if (0)
			param.solver = @ode23t;
			rk = Rietkerk(param);
			y0 = rk.init();
			tic
			[t,y] = rk.solve(y0,t);
			toc
			length(t)
	
			hold on
			plot(y(end,:))
end
if (0)
			% note this is
			opt = odeset('Jacobian',@rk.jacobian);
			t = linspace(0,param.T,ceil(param.T/0.5));
			tic()
			[t,y] = trapezoidal_fixed(@rk.dz_dt,t,cvec(y(10,:)),opt);
			toc()
			hold on
			plot(y(end,:))
end
