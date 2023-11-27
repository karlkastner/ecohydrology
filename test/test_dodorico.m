% Wed 18 Oct 12:29:57 CEST 2023
% Karl Kästner, Berlin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% note :
% a) without diffusion and constant coefficients, the value at each point converges to either value
%        of the two homogeneous stable states
% b) which state is chosen, is decided by the initial value, separated by the pole
%	hs1 ... pole  ... hs2
% c) with constant initial value, patterns can only form when the pole spatially varies
%    the poles are only influenced by l0 and b
% d) with diffusion, dx -> 0 and sd/dx const (gaussian error propagation)
%    the stationary points move towards the pole, until the distribution
%    of point values becomes uniform, i.e. there is no more pattern (bimodal distribution)
% -> since discretizing is equal to low-pass filtering, patterns only form when
%    the randomly varied parameter or ic is low-pass filtered
% e) dodorico 2006 and 2007 have different definitions of f

% implement dto !!!
% implement output class
% implement splitting scheme
%

T  = 1000;
dt = 1;
dt = 0.5;
%k = 0.5;
dx = 1;
L  = 100*[1,1];

nx  = L./dx;
%nx = 200*[1,1];
%L  = 0.25*[100,100];
mf  = 0.00;
sf  = 0.0;
dto = 0.5;
e   = 0.3;
rad = Dodorico('T',T,'L',L,'nx',nx,'pmu.f',mf,'pss.f',sf,'opt.path_str','mat/','opt.dto',dto ...
		,'opt.dt', dt ...
		,'pmu.ex',e ...
		,'pmu.ey',e ...
		,'pmu.w0',0.04 ...
		,'opt.solver', 'euler-forward' ...
		,'opt.mode', 'discrete' ...
		);
zp = rad.poles();
%rad.initial_condition = struct('mu',zp,'sd',sqrt(eps),'dist','uniform');
zp = 1.0;
rad.initial_condition = struct('mu',zp+sqrt(eps),'sd',0,'dist','uniform');
[t,z] = rad.run();
v = rad.extract2(z);

clf
subplot(2,2,1)
[x,y] = rad.x;
imagesc(x,y,squeeze(v(end,:,:)));
colorbar
axis equal
axis tight
subplot(2,2,2)
zbar = mean(z,2);
zsd  = std(z,[],2);
d = (rms(diff(z),2));
d(:,2) = max(abs(diff(z)),[],2);
%d = inner2outer(d);
plot(mid(t),d)
%plotyy(t,[zbar,zsd],t,d)
yyaxis right
plot(t,[zbar,zsd])

subplot(2,2,3)
hist(flat(v(end,:,:)))
vline(mean(v(end,:,:),'all'))
v_ = squeeze(v(end,:,:));




Shat = abs(fft2(v_-mean(v_,'all'))).^2;

[Sr,fr] = periodogram_radial(Shat,rad.L);
subplot(2,2,4);
plot(fr,Sr.normalized,'.-');

