function [delu] = dmccalc(Sp,Kmat,sn,delup,d,r,u,k,n)
%
% for use with dmcsim.m
% b.w. bequette
% 2 oct 00
% calculate the optimum control move
%
% first, calculate uold = u(k-n+1)...u(k-n+p)
%
  [m,p] = size(Kmat);
    uold = zeros(p,1);
  for i = 1:p;
       if k-n+i>0;
       uold(i) = u(k-n+i);
       else
       uold(i) = 0;
       end
  end
  dvec  = d*ones(p,1); % distúrbio aditivo
  rvec   = r*ones(p,1); % vetor de set-point
  y_free = Sp*delup + sn*uold + dvec;
  e_free = rvec-y_free;
  delu = Kmat(1,:)*e_free; % Obs.: Kmat(1,:), pois somente a primeira linha da matriz Kmat é implmentada