
clear all

%% Input data
% Set up model and discretize:
s = tf('s');
g11 = 1/(s*(s+1))*10/(10*s+1);
g21 = 1/(s*(s+10));
g31 = 1/(s*(s^2+4*s+5));
g12 = 1/s;
g22 = 1/(s*(s+7));
g32 = 5/(s*(s+5));
G = [g11 g12; g21 g22; g31 g32].';
MPC_case.Ts = 0.2; % sampling time
MPC_case.TF = c2d(G,MPC_case.Ts);
[MPC_case.A, MPC_case.B, MPC_case.C, MPC_case.D] = ssdata(MPC_case.TF);
% MPC_case.A = [1 3 2; 0 1 2; 0 0 1];
% MPC_case.B = [1 0; 2 0; 1 1];
% MPC_case.C = [1 0 0; 0 1 0];
% MPC_case.D = [0 0; 0 0];

MPC_case.nx = length(MPC_case.A);
[MPC_case.ny,MPC_case.nu] = size(MPC_case.D);
MPC_case.C_pinv = pinv(MPC_case.C);

% Weightings:
MPC_case.Qy = eye(MPC_case.ny);  % Output weight -> state weight for 1<i<Np-1;
MPC_case.Sy = 10*eye(MPC_case.ny);  % Terminal weight at i=Np
MPC_case.Qu = 1*eye(MPC_case.nu); % Input weight for 0<i<Np-1

% Constraints (implemented as hard here):
MPC_case.umin = -40*ones(MPC_case.nu,1); 
MPC_case.umax = 40*ones(MPC_case.nu,1);
MPC_case.ymin = -50*ones(MPC_case.ny,1); 
MPC_case.ymax = [50; 50];
MPC_case.delta_u_min = -5*ones(MPC_case.nu,1); 
MPC_case.delta_u_max = 20*ones(MPC_case.nu,1);
MPC_case.xmin = [-1e6*ones(MPC_case.nx,1); MPC_case.umin];
MPC_case.xmax = [1e6*ones(MPC_case.nx,1); MPC_case.umax];

% References for outputs, states and inputs:
MPC_case.yref_all = @(t)[(t<20)*40+(t>=20)*10; (t<25)*10+(t>=25)*49];     % Reference for 0<i<Np-1

% Horizon is Np (both for objective function and for constraints);
MPC_case.Np = 15;

% assign x0 and t_end
MPC_case.x0 = zeros(MPC_case.nx+MPC_case.nu,1); % First state value is zero
MPC_case.t_end = 40;

% Noise
MPC_case.noise_order = 100;
MPC_case.noise_percent = 0;

%% Run MPC
[U,X,Y,Y_measured] = MPC_simulation_delta_u(MPC_case);

%% plot results
t = 0:MPC_case.Ts:MPC_case.t_end;
t = [-1 t];
y1 = [0 0 Y(1,:)];
y2 = [0 0 Y(2,:)];
yref = MPC_case.yref_all(t);
yref(:,1) = [0; 0];

figure(537); clf
subplot(511),
plot(t,y1,[-1 MPC_case.t_end],[MPC_case.ymax(1) MPC_case.ymax(1)],'g--',[-1 MPC_case.t_end],[MPC_case.ymin(1) MPC_case.ymin(1)],'g--')
hold on;
stairs(t,yref(1,:),'r');
ylabel('Y_1')
axis([-1, MPC_case.t_end, min(y1)-1, max(y1)+1])

subplot(512),
plot(t,y2,[-1 MPC_case.t_end],[MPC_case.ymax(2) MPC_case.ymax(2)],'g--',[-1 MPC_case.t_end],[MPC_case.ymin(2) MPC_case.ymin(2)],'g--')
hold on;
stairs(t,yref(2,:),'r');
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

subplot(515)
U3 = [0 U(3,:) U(3,end)];
stairs(t,U3)
axis([-1 MPC_case.t_end MPC_case.umin(3)-0.3 MPC_case.umax(3)+0.3])
ylabel('u_3')

xlabel('time [s]')



