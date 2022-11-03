% Sat 11 Dec 16:10:32 CET 2021
	load('mat/rietkerk-1565959574.mat','t','y','rk');
	
	x = rk.x;
	b = rk.extract1(y);

	figure(1);
	clf
	subplot(2,2,1)
	imagesc(rk.x,t,b);

	subplot(2,2,2)
	plot(rk.x,b(end-1:end,:))

	c = rk.celerity(double(y));
	c_ = migration_celerity_1d(b,t(2)-t(1),x(2)-x(1));
	subplot(2,2,3)
	cla
	plot(t,c);
	hold on
	plot(mid(t),c_)
	ylim([-0.05,0.00])

	subplot(2,2,4)
	n=1;
	d=0.01;
	dt= t(2)-t(1)
	x_ = rk.dx*(-n:d:n);
	y_ = interp1(-n:n,xcorr(b(end,:),b(end-1,:),n),-n:d:n,'spline');
	plot(x_/dt,y_/max(y_));
	[mv,mdx]=max(y_);
	vline(x_(mdx)/dt)
	x_(mdx)/dt
	vline(-c(end),'linestyle','--','color','r')
	cmu = mean(c(end-100:end)), rk.L/cmu, 100/cmu

	
