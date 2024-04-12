% 2023-12-08 13:52:01.016536346 +0100

syms b w;
z = [b; w];
rad = Klausmeier();
rad.make_symbolic();
dz_dt = rad.dz_dt_react(0,z);
syms J;
J = [diff(dz_dt,z(1)),diff(dz_dt,z(2))]
J_ = rad.jacobian_react(0,z)
% test
J-J_

