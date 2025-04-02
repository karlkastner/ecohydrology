% 2024-12-17 13:16:08.553009225 +0100
% Karl Kastner, Berlin
dx  = 1;
nx  = 128*[1,1];
nn = prod(nx);
nvar =2;
L   = nx*dx;
rad = RAD_Linear();
rad.nvar = nvar;
rad.nx = nx;
rad.L  = L;
d = [0.1, 100];
c = [0.05, 1; % 0.1
     -1, 50]; % 50
rad.p.ex = d;
rad.p.ey = d;
rad.p.vx = 0.*d;
rad.p.vy = 0.*d;
rad.init_advection_diffusion_matrix();

[fx,fy,frr] =fourier_axis_2d(L,nx);

% direct perturbation of z (simply take impulse here)
%e = zeros(nn*nvar,1);
%e = e + 0.001*randn(nn*nvar,1);
%e(1) = 1;
e0 = normpdf(frr,0,0.02);
e0 = flat(e0);
e = [e0;e0];
%e = e0+circshift(e0,[30,0]);
%e = e+circshift(e0,[0,20]);
%5e = e+circshift(e0,[10,40]);
%e = [flat(e); zeros(nn,1)];
%e = e+[flat(normpdf(frr,0.2,0.02)); zeros(nn,1)];
%e = [flat(normpdf(frr,0,0.1)); zeros(nn,1)];
e = e + 0.1*randn(2*nn,1);
D2 = rad.aux.D2x + rad.aux.D2y;
I = speye(nn);
A = [c(1,1)*I + d(1)*D2, c(1,2)*I;
     c(2,1)*I, c(2,2)*I + d(2)*D2];

z = -A \ e;
z1 = reshape(z(1:nn),nx);

[cfit,A] = rad.fit(z,e,[],true);
cfit
[cfit,D2f] = rad.fit_fourier(z,e);
cfit

subplot(2,2,1);
imagesc(fftshift(z1))
axis square
subplot(2,2,2)
plot(fftshift(z1(1,:)))
subplot(2,2,3)
d2z1  = reshape(D2*z1(:),nx);
%imagesc(d2z1)
plot(d2z1(1,:))
subplot(2,2,4)
d2z1_ = real(ifft2(reshape(D2f,nx).*fft2(z1)));
%imagesc(d2z1_)
plot([d2z1_(1,:)'-d2z1(1,:)',d2z1(1,:)'])


