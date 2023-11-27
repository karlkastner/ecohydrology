% 2023-10-16 16:09:01.485986224 +0200
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

		T   = 40;
		dto = 10;
		dx = 1;
		L  = 1*[100,100];
		nx = L/dx;
		dt = 0.01;
		vw = [0,0];
		ew = [300,300];


		 k = Klausmeier('nx',nx,'T',T,'L',L,'opt.solver','solve_split','opt.dt',dt ...
			,'pmu.vw',vw ...
			,'pmu.ew',ew ...
		        ,'opt.dto', dto ...
			...'pmu.eb',eb, ...
			...'psd.g', s, ...
			... % 'pmu.r', 0.115, ...
			... % 'pmu.r', 0.15, ...
			...'pmu.r', 0.96 ...
			);
%, ...			 'odeopt',oo);
		 k.init();
		 tic();
		 [t,zz] = k.solve();
		 toc
		 b1 = k.extract2(zz(:,end));

		clf
		 subplot(2,3,1);
		 imagesc(b1);
		colorbar

		 oo=odeset('reltol',1e-4,'abstol',1e-4);
		 k = Klausmeier('nx',nx,'T',T,'L',L,'opt.solver',@ode23,'opt.dt',dt ...
			,'pmu.vw',vw ...
			,'pmu.ew',ew ...
		        ,'opt.dto', dto ...
			...'pmu.eb',eb, ...
			...'psd.g', s, ...
			... % 'pmu.r', 0.115, ...
			... % 'pmu.r', 0.15, ...
			...'pmu.r', 0.96 ...
			 ,'odeopt',oo ...
			);
		 k.init();
		 tic();
		 [t,zz] = k.solve();
%		[x,y] = k.x;
%		x=1;
%		y=1;
		toc()
		%hold on
		b2 = k.extract2(zz(end,:));

%		z = zz(end,1:prod(nx));
%		z = reshape(z,k.nx);
%		clf();
		 subplot(2,3,2);
		 imagesc(b2);
		colorbar

		subplot(2,3,3)
		imagesc((b2-b1));
		colorbar
		rms(flat(b2-b1),'all')./rms(b2,'all')

if (0)
		S = abs(fft2(z-mean(z,'all'))).^2;
		subplot(2,3,2);
		imagesc(fftshift(S))		
		
		[Sr,fr] = periodogram_radial(S,k.L);
		subplot(2,3,3);
		plot(fr,Sr.normalized,'.-');
		
		
		subplot(2,3,4)
		bbar = mean(zz(:,1:end/2),2);
		plot(t,bbar)
end
		 %plot(z(end,1:nx));  
