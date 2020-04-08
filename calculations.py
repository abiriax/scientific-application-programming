rho_0 = 1.64e-24
p = 3/5
v_0 = 1e5
p_cgs = p * v_0**2 * rho_0
l = 1.49597871e+13
M_0 = rho_0*l**3
M_sun_cgs = 2.e33
M_code = (6.67e-8)**-1
M_cgs = M_0 * M_code 
M_dot = 1.4e11*(M_cgs/M_sun_cgs)**2 *(rho_0/1.e-24)*(v_0/1.e6)**-3

print(M_dot)

