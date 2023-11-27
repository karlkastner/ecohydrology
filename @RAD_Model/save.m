% Mon 31 May 20:20:46 CEST 2021
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
function save(rad,t,y,out)
	% clear auxiliary variables to save disk space
	rad.aux = struct();
	[oname,oname_final] = rad.filename();
	% write model results to hard drive
	save(oname,'-v7.3','t','y','rad','out');
	if (length(t)>1)
		t = t([1,end]);
		y = [cvec(rad.z0),cvec(out.y_final)];
	end
	save(oname_final,'-v7.3','t','y','rad','out');
end

