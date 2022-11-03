% Tue  3 May 13:12:57 CEST 2022

			% 2d
			param = struct();
			param.L = 100*[1,1];
			n       = 100;
			param.n = [n,n+2];
			dt = 0.01;
			param.T = [10];
			%param.nt = 2;
			param.dto = 1;
			param.dt = dt;
			param.pss.a = 0; % 0.1;

			param.pmu.eh = [0,100];
			param.pmu.vh = [10,0];

			rk = Rietkerk(param);
			y0 = rk.init();
			nf = 5;

if (1)
			[b,w,h]=rk.extract2(rvec(y0));
			b = squeeze(b);
			w = squeeze(w);
			h = squeeze(h);
			b = trifilt1(b,nf);
			w = trifilt1(w,nf);
			h = trifilt1(h,nf);
			y0 = [b(:);w(:);h(:)];
end
			tic
%			rk.odeopt = odeset('MaxStep',1e-3);

if (~exist('b_ref','var') || t(end) ~= param.T)
			[t,y] = rk.solve(y0);
			dt_ref = t(end)/length(t)
			tr_ref = toc
			b = rk.extract2(y);
			b1 = squeeze(b(1,:,:));
			b = squeeze(b(end,:,:));
			b_ref = b;
end
			figure(1);
			subplot(2,2,1);
			cla
			%plot(b1(end/2,:))
			%hold on
			plot(b_ref(end/2,:));
			hold on
			subplot(2,2,2)
			cla
			plot(b_ref(:,end/2));
			hold on

			%imagesc(
if (0)			
			[t,y] = euler(@rk.dz_dt,t,y0);
			b = rk.extract2(y);
			b = squeeze(b(end,:,:));
			subplot(2,2,1);
			plot(b(end/2,:))
			subplot(2,2,2)
			plot(b(:,end/2))
end
	
			%dt = [0.01 0.02 0.05 0.1 0.2,0.5,1];
			%dt = [0.05 0.1 0.2,0.5,1,2];
			dt = [0.1,0.2,0.5,1,2];
			%dt  = 1;
			%dt = [0.1] % 0.02 0.05 0.1 0.2,0.5,1,2];
			tr = [];
			res = [];
			res1 = [];
			%method = [0,1,2];
			method = 2;
			for jdx=1:length(method)
			for idx=1:length(dt)	
			fallback = [];	
				t = linspace(0,param.T,round(param.T/dt(idx))+1);
					tic();
					[t,y,fallback(:,jdx)] = rk.solve_split(t,y0,dt(idx),method(jdx));
					%[t,y,fallback(:,jdx)] = rk.solve_split(y0,t,method(jdx));
					tr(idx,jdx) = toc();
					bb = rk.extract2(y);
					b(:,:,jdx) = squeeze(bb(end,:,:));
					b(:,:,idx) = squeeze(bb(end,:,:));

			if (1 == idx)
				figure(1);
				subplot(2,2,1);
				if (1==jdx) cla; end
				plot(b(end/2,:,jdx));
				subplot(2,2,2)
				plot(b(:,end/2,jdx));
				%bx(:,idx) = b(:,end/2);
			end
				if (idx == 1)
					b1(:,:,jdx)= b(:,:,jdx);
				end
				res(idx,jdx) = rms(flat(b(:,:,jdx)-b_ref));
				res1(idx,jdx) = rms(flat(b(:,:,jdx)-b1(:,:,jdx)));
			if (1==jdx)
				figure(2)
				if (1==idx) clf; plot(b_ref(1,:)); end
				hold on
				plot(squeeze(bb(end,1,:)))
			end
			end
			end
			res = res./rms(flat(b_ref));
			res1 = res1./rms(flat(b_ref));
			disp(tr_ref)
			disp(tr)

			%subplot(2,2,3)
			%b = squeeze(bb(1:10:end,end/2,:));
			%plot(squeeze(b)')

			figure(1)
			subplot(2,3,4);
			cla
			plot(t,cumsum(fallback));
	

			if (length(dt)>1)
%			c=[ones(length(dt)-1,1),log(cvec(dt(2:end)))] \ log(res(2:end,2)); hold on;

			subplot(2,3,5)
			cla
			loglog(dt,res,'.-')
			hold on
			loglog(dt,res1,'.--')
%			plot(dt,exp(c(1)).*dt.^c(2))
%			yyaxis right
			subplot(2,3,6)
			cla
			loglog(dt,tr,'.-')
			end
			

			
			
			
		
			


