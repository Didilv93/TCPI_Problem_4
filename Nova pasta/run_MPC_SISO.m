
clear all

%% Input data
% Set up model and discretize:
s = tf('s');
G = 1/(s*(s+1))*10/(10*s+1);
MPC_case.Ts = 0.2; % sampling time
MPC_case.TF = c2d(G,MPC_case.Ts);
[MPC_case.A, MPC_case.B, MPC_case.C, MPC_case.D] = ssdata(MPC_case.TF);
MPC_case.nx = length(MPC_case.A);
[MPC_case.ny,MPC_case.nu] = size(MPC_case.D);
MPC_case.C_pinv = pinv(MPC_case.C);

% Weightings:
MPC_case.Qy = 1;  % Output weight -> state weight for 1<i<Np-1;
MPC_case.Sy = 10;  % Terminal weight at i=Np
MPC_case.Qu = .01; % Input weight for 0<i<Np-1

% Constraints (implemented as hard here):
MPC_case.umin = -30; 
MPC_case.umax = 30;
MPC_case.ymin = -50; 
MPC_case.ymax = 50;

% References for outputs, states and inputs:
MPC_case.yref = 10;     % Reference for 0<i<Np-1
MPC_case.uref = 0;     % Input reference is here assumed zero

% Horizon is Np (both for objective function and for constraints);
MPC_case.Np = 20;

% assign x0 and t_end
[A,B,C,D] = ssdata(MPC_case.TF);
nx = length(A);
MPC_case.x0 = zeros(nx,1); % First state value is zero
MPC_case.t_end = 25;

% Noise
MPC_case.noise_order = 100;
MPC_case.noise_percent = 0;

%% Run MPC
[U,X,Y,Y_measured,MPC_case] = MPC_simulation(MPC_case);

%% plot results
t = 0:MPC_case.Ts:MPC_case.t_end;
t = [-1 t];
y = [0 0 Y];
y_meas = [0 0 Y_measured];

figure(532)
subplot(311),
plot(t,y,[-1 0],[0 0],'r:',[0 0],[0 MPC_case.yref],'r:',[0 MPC_case.t_end],[MPC_case.yref MPC_case.yref],...
      'r:',[-1 MPC_case.t_end],[MPC_case.ymax MPC_case.ymax],'g--',[-1 MPC_case.t_end],[MPC_case.ymin MPC_case.ymin],'g--')
ylabel('output y')
axis([-1, MPC_case.t_end, min(y)-0.2, max(y)+0.2])

subplot(312),
plot(t,y_meas,[-1 0],[0 0],'r:',[0 0],[0 MPC_case.yref],'r:',[0 MPC_case.t_end],[MPC_case.yref MPC_case.yref],...
      'r:',[-1 MPC_case.t_end],[MPC_case.ymax MPC_case.ymax],'g--',[-1 MPC_case.t_end],[MPC_case.ymin MPC_case.ymin],'g--')
ylabel('Measured y')
axis([-1, MPC_case.t_end, min(y_meas)-0.2, max(y_meas)+0.2])

subplot(313)
U = [0 U U(end)];
stairs(t,U)
axis([-1 MPC_case.t_end MPC_case.umin-0.1 MPC_case.umax+0.1])
ylabel('input u')
xlabel('time [s]')

