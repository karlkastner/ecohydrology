% 2024-12-22 10:04:26.779823832 +0100
% Karl Kastner, Berlin

nvar = 1;
n = 128*[1,1];
L = n;
d=[1.3,1.4,1.5];
d=[2,3,4];
%d=[0,0,0];
d = d(1:nvar);
rng(1);

x = innerspace(0,1,n(1));
%x = innerspace(0,1,n);
r = sin(2*pi*x);
r = cos(2*pi*x);

A = eye(nvar)+rand(nvar);
%A = eye(nvar);
%A = diag([1,10]);
%A(1,2) = -2
%A = 1./(2:10).^2;
A =reshape(A(1:nvar.^2),nvar,nvar);
if (2 == length(n))
rhs1 = flat(cvec(r)*rvec(r));
else
	rhs1 = cvec(r);
end

%rhs=randn(n*n,nvar);
rhs        = rhs1;
if (nvar>1)
rhs(:,2) = rhs1;
end
if (nvar>2)
rhs(:,3) = rhs1;
end
%.^3;
%rhs(:,1) = 0; 
%rhs(:,2) = 0; 
%rhs(:,3) = 0; 
rhs = rhs(:);

if (2 == length(n))
I = speye(n(1)*n(2));
[Dx,Dy,Dxx,Dxy,Dyy] = derivative_matrix_2d(n,L./n,2,{'circular','circular'},true);
D2 = Dxx+Dyy;
else
	I = speye(n);
	D2 = derivative_matrix_2_1d(n,L/n,2,'circular','circular',true);
end

AA = kron(A,I) + kron(diag(d(1:nvar)),D2);
z = AA\rhs(:);
%dz_dt = A*z;

%flat(ones(n*n,1)*[2,3,5]);
%rhs=reshape(rhs,n,n,nvar);
%A=A';
z_  = squeeze(rad3_solve_fourier(n,L,A,d,rhs,nvar));

disp('testing z-z');
rms(z(:)-z_(:))
rms(reshape(z(:)-z_(:),[],nvar))

dz_dt = radlin_dz_dt_fourier(n,L,A,d,z,nvar);
dz_dt_ = radlin_dz_dt_fourier(n,L,A,d,z_,nvar);
disp('testing ||A z - rhs|| = 0')
rms(reshape(AA*z(:)-rhs(:),[],nvar))
rms(reshape(AA*z_(:)-rhs(:),[],nvar))
rms(reshape(dz_dt(:)-rhs(:),[],nvar))
rms(reshape(dz_dt_(:)-rhs(:),[],nvar))

%z__= squeeze(rad3_solve_fourier([n,n],[1,1],A,d,rhs));

%rms(z(:)-z__(:))
%rms(z_(:)-z__(:))

if (length(n) == 1)
rhs  = reshape(rhs,[],nvar);
z  = reshape(z,[],nvar);
z_ = reshape(z_,[],nvar);
clf
for idx=1:nvar
	subplot(2,nvar,idx)
	plot([rhs(:,idx),z(:,idx),z_(:,idx)]);
end
end
