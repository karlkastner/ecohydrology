% Wed 18 Oct 08:42:39 CEST 2023
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
% fisher 1939, Murray, 1989 for determinisitic part
function dz_dt = dz_dt_react(obj,t,z)
	% note 1 : F in Dodorico 20007 models the death by fire as a Poisson
	% process with rate l = l0 + b*V, i.e. on average occurs a fire
	% event in a particular grid cell every l time steps,
	% here, a continuous approach is chosen,
	% where at every time step vegetation proportional to 1/l dies off 
	%	  note that one of the stable states is larger zmax,
	%	  which should be interpreted as a retardation value rather than maximum
	% note 2 : dodorico 2007 does not mention a time or space step,
	%          though his choice of exponential and poisson distribution
	%	   seem to imply a time step of 1, and a space step of 1
	%          however,  system is numerically unstable for the parameters
	%	   proposed by DOrico:
	%	   dt_max <= dx^2/(4*e) = 1^2/(4*0.3) = 0.83
	% note 3: with the poisson distribution, several events can occur in the same time step,
	%	  as the pdf of the poissin distribution larger 0 at the origin
	%	  dodorico does not specify how this is treated

	% note that dodorico substracts f (inverse definition of sign)
	l = obj.p.l0 + obj.p.b.*z;
	%l = obj.p.l0.*exp(obj.p.b.*z);

	switch (obj.opt.mode)
	case {'continuous'}
		F  = obj.p.w0./l;
	case {'discrete'}
		% note that dodorico chooses the iverse definition as matlab for
		% the poisson distribution parameter
		[event,obj.aux.next] = poisson_noise(t/obj.opt.dt,obj.aux.next,1./l);
		%[median(l),median(event)]
		mu = obj.p.w0/obj.opt.dt;
		sd = obj.p.w0/sqrt(obj.opt.dt);
		[a,b] = gamma_moment2par(mu,sd);
		%F = exprnd(obj.p.w0,prod(obj.nx),1);
		F = gamrnd(a,b).*event;
	end
	if (~issym(z))
		dz_dt = obj.p.a.*(max(z,0)+obj.p.ze).*(obj.p.zmax-z) - F.*(z>0);
	else
		dz_dt = obj.p.a.*(z+obj.p.ze).*(obj.p.zmax-z) - F;
	end
end

