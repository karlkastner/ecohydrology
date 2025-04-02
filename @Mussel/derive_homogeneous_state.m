rk.make_symbolic; syms m a; dz_dt = rk.dz_dt_react(0,[m,a]); a_ = solve(dz_dt(1),a), m_ = solve(subs(dz_dt(2),a,a_),m), a_=solve(subs(dz_dt(1),m,m_),a)
