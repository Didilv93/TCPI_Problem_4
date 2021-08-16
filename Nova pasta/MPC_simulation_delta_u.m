
% This function simulates MPC with iterating through time, and each time
% implementing the first (in time) input vector that is solved by 
% optimization from MPC_calculation.m
%
% The formulation and original code (for SISO systems) is by Elling W.
% Jacobsen from KTH university, Sweden. The formulation can be found here:
% http://www.kth.se/polopoly_fs/1.120161!/Menu/general/column-content/attachment/f13_ewj.pdf
% The code is modified and generalized for MIMO systems by Pooya Rezaei, 
% University of Vermont, USA. 

function [U,X,Y,Y_measured,MPC_case] = MPC_simulation_delta_u(MPC_case)

Nsim = MPC_case.t_end/MPC_case.Ts; % The number of simulation intervals
% Initialize variables
x0 = MPC_case.x0;
x_noisy = MPC_case.x0;
U = zeros(MPC_case.nu,Nsim);
X = zeros(MPC_case.nx,Nsim);
Y = zeros(MPC_case.ny,Nsim);
Y_measured = zeros(MPC_case.ny,Nsim);

for k = 1:Nsim
    t_Np = (k-1)*MPC_case.Ts:MPC_case.Ts:(k-1+MPC_case.Np-1)*MPC_case.Ts;
    % References for outputs and inputs:
    yref_Np = MPC_case.yref_all(t_Np);
    % Calculate future Np inputs
    [u,MPC_case] = MPC_calculation_delta_u(x_noisy,yref_Np,MPC_case);
    % Implement present input for the plant model
    x_out = MPC_plant(x0(1:MPC_case.nx),u,MPC_case);
    x_noisy = Addnoise(x_out,MPC_case.noise_order,MPC_case.noise_percent);
    % Save states and set initial value for next time slot
    U(:,k) = u;
    Y(:,k) = MPC_case.C*x_out;
    Y_measured(:,k) = MPC_case.C*x_noisy;
    X(:,k) = x0(1:MPC_case.nx);
    x0 = [x_out;u];
    x_noisy = [x_noisy;u];
end

