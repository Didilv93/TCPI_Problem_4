function x_out = MPC_plant(x0,u,MPC_case)

x_out = MPC_case.A*x0+MPC_case.B*u; 
