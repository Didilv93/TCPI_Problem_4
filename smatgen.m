function [Sf,Sp,Kmat] = smatgen(s,p,m,n,w)
%
% b.w. bequette
% 28 Sept 00, revised 2 Oct 00
% generates dynamic matrix and feedback gain matrix
% assumes s = step response column vector
% Sf   = Dynamic Matrix for future control moves (forced)
% Sp   = Matrix for past control moves (free)
% Kmat = DMC feedback gain matrix
% s    = step response coefficient vector
% p    = prediction horizon
% m    = control horizon
% n    = model horizon
% w    = weight on control input
%
% first, find the dynamic matrix
  for j = 1:m
         Sf(:,j) = [zeros(j-1,1);s(1:p-j+1)];       
  end
%
% now, find the matrix for past moves
%
  for i = 1:p
         Sp(i,:) = [s(i+1:n-1)' zeros(1,i-1)];
  end
%
% find the feedback gain matrix, Kmat
%
  Kmat = inv(Sf'*Sf + w*eye(m))*Sf';
  
  
  