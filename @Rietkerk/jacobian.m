% 2021-07-05 16:48:44.611627340 +0200
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
%% jacobian of the Rietkerk model
%
function [A,res] = jacobian(obj,t,z,withda)
	if (nargin()<4)
		withda = true;
	end

	p = obj.p;
	[b,w,h] = obj.extract1(z);

	U_div_b  = p.gb.*w./(w + p.kw);
        In = p.a.*h.*(b + p.kb.*p.w0)./(b+p.kb);

	if (withda)
		A = [ p.eb*obj.aux.D2x + (p.cb.*U_div_b - p.db).*obj.aux.I, ((b.*p.cb*p.gb*p.kw)./(p.kw + w).^2).*obj.aux.I,     obj.aux.Z
			((p.a.*h)./(b + p.kb) - (p.gb.*w)./(p.kw + w) - (p.a.*h.*(b + p.kb.*p.w0))./(b + p.kb).^2).*obj.aux.I, ...
				p.ew*obj.aux.D2x + (-p.rw  - (b.*p.gb)./(p.kw + w) + (b.*p.gb.*w)/(p.kw + w).^2).*obj.aux.I, ...
				(In./h).*obj.aux.I;
			(p.a.*h.*p.kb.*(p.w0 - 1))./(b + p.kb).^2.*obj.aux.I,    obj.aux.Z, p.eh*obj.aux.D2x + p.vh*obj.aux.D1x - (In./h).*obj.aux.I
		];
	else
		% reaction part only
	A = [ + (p.cb.*U_div_b - p.db).*obj.aux.I, ((b.*p.cb*p.gb*p.kw)./(p.kw + w).^2).*obj.aux.I,     obj.aux.Z
		((p.a.*h)./(b + p.kb) - (p.gb.*w)./(p.kw + w) - (p.a.*h.*(b + p.kb.*p.w0))./(b + p.kb).^2).*obj.aux.I, ...
			(-p.rw  - (b.*p.gb)./(p.kw + w) + (b.*p.gb.*w)/(p.kw + w).^2).*obj.aux.I, ...
			(In./h).*obj.aux.I;
		(p.a.*h.*p.kb.*(p.w0 - 1))./(b + p.kb).^2.*obj.aux.I,    obj.aux.Z, - (In./h).*obj.aux.I
	];
	end
	if (nargout()>1)
	res = [ (b.*p.cb*p.gb.*h.*p.kb)./(h + p.kb).^2
	       b.*( -p.gb.*h.*p.kb./(h + p.kb).^2 + p.a.*w.*p.kw.*(1 - p.w0)./(b + p.kw).^2 )
              -p.R - p.a.*b.*p.kw.*w.*(1 - p.w0)./(b + p.kw).^2];
	end
end

