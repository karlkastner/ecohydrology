% Fri  8 Dec 13:50:59 CET 2023
function J = dz_dt_react(obj,t,z)
	[b,w] = obj.extract1(z);
	c = obj.p.c;
	g = obj.p.g;
	l = obj.p.l;
	d = obj.p.d;
	if (issym(z))
	J = [2*b*c*g*w - d,     b^2*c*g
	          -2*b*g*w, - g*b^2 - l];
	else
	J = [diag(sparse(2*c*g*b.*w - d)),     diag(sparse(c*g*(b.*b)))
	          diag(sparse(-2*g*(b.*w))), diag(sparse(- g*(b.*b) - l))];
	end
end

