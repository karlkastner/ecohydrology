 syms a m;
 rk.make_symbolic();
 dz_dt = rk.dz_dt_react(0,[m,a]);
 J = [diff(dz_dt,m),diff(dz_dt,a)]
