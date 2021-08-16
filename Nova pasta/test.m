
clear all

%% Input data
% Set up model and discretize:
s=tf('s');
Ts = 0.5;
N = 50;
Np = 10;
Rh = 6.600000e+00;
Ch = 2.000000e-01;
CT = 3.000000e+00;
RT = 3.250000e-01;
A  = 2.500000e-02;
h_ = 3.000000e-01;
T_ = 3.300000e+01;
Ti_ = 1.800000e+01;
T2_ = 3.600000e+01;
Fc_ = 2.000000e-09;
Fi_ = 1.000000e-02;
qi_ = 1.200000e-08;

a = -((Fi_ + Fc_)/A*h_ - 1/(CT* RT));
b = -((Fi_*(Ti_ - T_) + Fc_*(T2_ - T_))/(A*(h_^2)));
c = 1/CT;
d = (T2_ - T_)/A*h_;
e = (Ti_ - T_)/A*h_;
f = (Fi_/(A*h_) - 1/(CT*RT));
g = Fc_/(A*h_);

A = [(-1)/(Rh*Ch) 0; b a];
B = [1/Ch 0; d c];
C = [0 0; b/(s + a) 0];
D = [Rh/(Rh*Ch*s + 1) 0; d/(s + a) c/(s + a)];

G_s = (C*(s*eye(2) - A)^-1)*B + D;

MPC_case.Ts = Ts; % sampling time
MPC_case.TF = c2d(G_s,MPC_case.Ts);
[MPC_case.A, MPC_case.B, MPC_case.C, MPC_case.D] = ssdata(MPC_case.TF);

MPC_case.nx = length(MPC_case.A);
[MPC_case.ny,MPC_case.nu] = size(MPC_case.D);
MPC_case.C_pinv = pinv(MPC_case.C);

% Weightings:
MPC_case.Qy = eye(MPC_case.ny);  % Output weight -> state weight for 1<i<Np-1;
MPC_case.Sy = Np*eye(MPC_case.ny);  % Terminal weight at i=Np
MPC_case.Qu = .2*eye(MPC_case.nu); % Input weight for 0<i<Np-1

% Constraints (implemented as hard here):
MPC_case.umin = -40*ones(MPC_case.nu,1); 
MPC_case.umax = 40*ones(MPC_case.nu,1);
MPC_case.ymin = -50*ones(MPC_case.ny,1); 
MPC_case.ymax = [50; 18];

% References for outputs, states and inputs:
MPC_case.yref = [44; 14];     % Reference for 0<i<Np-1
MPC_case.uref = zeros(MPC_case.nu,1);     % Input reference is here assumed zero

% Horizon is Np (both for objective function and for constraints);
MPC_case.Np = 20;

% assign x0 and t_end
MPC_case.x0 = zeros(MPC_case.nx,1); % First state value is zero
MPC_case.t_end = 25;

% Noise
MPC_case.noise_order = 100;
MPC_case.noise_percent = 0;

%% Run MPC
[U,X,Y,Y_measured] = MPC_simulation(MPC_case);

%% plot results
t = 0:MPC_case.Ts:MPC_case.t_end;
t = [-1 t];
y1 = [0 0 Y(1,:)];
y2 = [0 0 Y(2,:)];

figure(534)
subplot(511),
plot(t,y1,[-1 0],[0 0],'r:',[0 0],[0 MPC_case.yref(1)],'r:',[0 MPC_case.t_end],[MPC_case.yref(1) MPC_case.yref(1)],...
      'r:',[-1 MPC_case.t_end],[MPC_case.ymax(1) MPC_case.ymax(1)],'g--',[-1 MPC_case.t_end],[MPC_case.ymin(1) MPC_case.ymin(1)],'g--')
ylabel('Y_1')
axis([-1, MPC_case.t_end, min(y1)-1, max(y1)+1])

subplot(512),
plot(t,y2,[-1 0],[0 0],'r:',[0 0],[0 MPC_case.yref(2)],'r:',[0 MPC_case.t_end],[MPC_case.yref(2) MPC_case.yref(2)],...
      'r:',[-1 MPC_case.t_end],[MPC_case.ymax(2) MPC_case.ymax(2)],'g--',[-1 MPC_case.t_end],[MPC_case.ymin(2) MPC_case.ymin(2)],'g--')
ylabel('Y_2')
axis([-1, MPC_case.t_end, min(y2)-1, max(y2)+1])

subplot(513)
U1 = [0 U(1,:) U(1,end)];
stairs(t,U1)
axis([-1 MPC_case.t_end MPC_case.umin(1)-0.3 MPC_case.umax(1)+0.3])
ylabel('u_1')

subplot(514)
U2 = [0 U(2,:) U(2,end)];
stairs(t,U2)
axis([-1 MPC_case.t_end MPC_case.umin(2)-0.3 MPC_case.umax(2)+0.3])
ylabel('u_2')

xlabel('time [s]')



