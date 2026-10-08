Sod-like shock tube along z (N=32, periodic, two interfaces at the
cell faces 0 and 16). Left (cells 1..16): rho = Dgas_avg = 1/6,
p = rho*c_L^2/gamma with c_L = 2.0; right: rho = 0.125*Dgas_avg, p = 0.1 p_L.
Uniform bulk velocity (vx,vy,vz) = (-0.15, 0.1, 0.25) (needed, see H1 in
docs/port_notes_hydro.md); the exact solution is the Riemann solution
advected with the bulk velocity. v_gas = 0 initially. Neutral w_gas.
Time step: the code's own Courant step (timestep = courant).
