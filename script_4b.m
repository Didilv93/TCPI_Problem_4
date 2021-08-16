s=tf('s');
Ts = 0.5;
N = 50;
tfinal = N * Ts;
t = 0:delt:tfinal;  % time vector
kfinal = length(t); % number of time intervals

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
G_z = c2d(G_s, Ts, 'zoh');

[phi, gamma, cd, dd] = ssdata(G_z);
[sys] = step(G_z,Ts:Ts:N*Ts);
p = 10;
m = 1;
weight = 0;
[Sf,Sp,Kmat] = smatgen(sys(:,1,1), p, m, N, weight);

% plant initial conditions
  xinit = zeros(size(A,1),1); % Ca-Cas = 0 e Cb-Cbs = 0
  uinit = 0; % valor inicial da variável manipulada
  yinit = 0; % valor inicial da resposta
% initialize input vector
  u     = ones(min(p,kfinal),1)*uinit;% vetor inicial u com número de linhas = min(p,kfinal)
%
  dup   = zeros(N-2,1); % vetor inicial delta u (ações passadas)
  sn    = sys(N,1,1);     % last step response coefficient
  x(:,1)= xinit;
  y(1)  = yinit;
  dist(1) = 0; % distúrbio aditivo
%
% set-up is done, start simulations
%
 for k = 1:kfinal
%
 du(k) = dmccalc(Sp,Kmat,sn,dup,dist(k),r(k),u,k,n);
% perform control calculation
  if k > 1
 u(k) = u(k-1)+du(k); % control input
  else
 u(k) = uinit + du(k);
  end
% plant equations
 x(:,k+1) = phi*x(:,k)+gamma*u(k);
 y(k+1) = cd*x(:,k+1);
% model prediction
  if k-n+1>0 
   ymod(k+1) = s(1)*du(k) + Sp(1,:)*dup + sn*u(k-n+1);
  else
   ymod(k+1) = s(1)*du(k) + Sp(1,:)*dup;
  end
% disturbance compensation
%
 dist(k+1) = y(k+1) - ymod(k+1);
% additive disturbance assumption
% put input change into vector of past control moves
 dup = [du(k);dup(1:n-3)];
 end
%
% stairs plotting for input (zero-order hold) and setpoint
%
  [tt,uu] = stairs(t,u);
  [ttr,rr] = stairs(t,r);
%
  figure(1)
  subplot(2,1,1)
  plot(ttr,rr,'--',t,y(1:length(t)))
  ylabel('y')
  xlabel('time')
  title('plant output')
  subplot(2,1,2)
  plot(tt,uu)
  ylabel('u')
  xlabel('time')
