% 2023-10-16 13:54:44.094332220 +0200
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


% n.b. celerity (phase speed) only correct for dt=0.01
if (1)
		T = 1000;
		nx = 250;
		dt = 0.1;

		ew = 300;
		vw = 0;

		

		oo=struct('reltol',1e-2,'abstol',1e-2);
		 k = Klausmeier('nx',nx,'T',T,'L',1e3,'opt.solver','solve_split','opt.dt',dt,'p.vw',vw,'p.ew',ew);
		 k.init();
		 tic();
		 [t,z] = k.solve();
		 toc
		 clf;
		 plot(z(1:nx,end));  
		 %k = Klausmeier('nx',nx,'T',T,'L',1e3,'opt.solver',@ode23,'opt.dt',0.1);
		 k = Klausmeier('nx',nx,'T',T,'L',1e3,'opt.solver',@ode23,'opt.dt',dt,'p.vw',vw,'p.ew',ew);
		 k.init();
		 tic();
		 [t,z] = k.solve();
		 toc
		 hold on
		 plot(z(end,1:nx));  
end
		

if (0)	
		 T = 1000;
		 nx = 250*4;
		 %dt = 0.01;
% ew : 300 : pattern
		 ew = 300;
		 ew = 100;
		 %sdb = 0.0045;
		 s = 1;
		 %s = 0;
		 k = Klausmeier('nx',nx,'T',T,'L',1e3,'opt.solver',@ode23,'opt.dt',0.1,'p.ew',ew,'p.vw',0,'ps.g',s);
%, 'p.r', 0.1);
		 k.init();
		 tic();
		 [t,z] = k.solve();
		 toc
		 clf
		subplot(2,2,1)
		 plot(z(end,1:nx));
		subplot(2,2,2) 
		S = cvec(abs(fft(z(end,1:nx)-mean(z(end,1:nx)))).^2);
		S(:,2)=abs(fft(k.p.g-mean(k.p.g))).^2;
		S = S./sum(S)
		plot(S)
end 
