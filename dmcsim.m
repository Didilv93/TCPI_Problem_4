% dmcsim
%
% 26 Oct 00 - revised 4/5 Aug 01 for Chapter 16 - MPC
% b.w. bequette
% unconstrained DMC simulation, calls: smatgen.m and dmccalc.m
% currently contains van de vusse reactor model and plant
%
% MPC tuning and simulation parameters
  n = 50;        % model length
  p = 10;        % prediction horizon
  m =  1;        % control horizon
  weight = 0.0;  % weighting factor
  ysp = 1;       % setpoint change (from 0)
  timesp = 1;    % time of setpoint change
  delt =   0.1;  % sample time
  tfinal = 6;    % final simulation time
  noise  = 0;    % noise added to response coefficients
%
  t = 0:delt:tfinal;  % time vector
  kfinal = length(t); % number of time intervals
  ksp = fix(timesp/delt);
  r = [zeros(ksp,1);ones(kfinal-ksp,1)*ysp]; % setpoint vector
%
% ----- insert continuous model here -----------
% model (continuous state space form)
%
  a = [-2.4048 0;0.8333 -2.2381]; % a matrix – van de vusse
  b = [7; -1.117];                % b matrix – van de vusse
  c = [0 1];                      % c matrix – van de vusse
  d = 0;                          % d matrix – van de vusse
  sysc_mod = ss(a,b,c,d); % create LTI "object"
%
% ----- insert plant here -----------
% perfect model assumption (plant = model)
  ap = a;
  bp = b;
  cp = c;
  dp = d;
  sysc_plant = ss(ap,bp,cp,dp);
%
% discretize the plant with a sample time, delt
%
  sysd_plant = c2d(sysc_plant,delt)
  [phi,gamma,cd,dd] = ssdata(sysd_plant)
%
% evaluate discrete model step response coefficients
%
  [s] = step(sysc_mod,[delt:delt:n*delt]);
%
% generate dynamic matrices (both past and future)
%
  [Sf,Sp,Kmat] = smatgen(s,p,m,n,weight);
%
% plant initial conditions
  xinit = zeros(size(a,1),1); % Ca-Cas = 0 e Cb-Cbs = 0
  uinit = 0; % valor inicial da variável manipulada
  yinit = 0; % valor inicial da resposta
% initialize input vector
  u     = ones(min(p,kfinal),1)*uinit;% vetor inicial u com número de linhas = min(p,kfinal)
%
  dup   = zeros(n-2,1); % vetor inicial delta u (ações passadas)
  sn    = s(n);     % last step response coefficient
  x(:,1)= xinit;
  y(1)  = yinit;
  dist(1) = 0; % distúrbio aditivo
%
% set-up is done, start simulations
%
 for k = 1:kfinal;
%
 du(k) = dmccalc(Sp,Kmat,sn,dup,dist(k),r(k),u,k,n);
% perform control calculation
  if k > 1;
 u(k) = u(k-1)+du(k); % control input
  else
 u(k) = uinit + du(k);
  end
% plant equations
 x(:,k+1) = phi*x(:,k)+gamma*u(k);
 y(k+1) = cd*x(:,k+1);
% model prediction
  if k-n+1>0;
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