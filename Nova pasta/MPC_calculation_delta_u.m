
% Sets up model and all matrices needed to construct Hessian and
% constraints in QP-problem. The matrices that do not change during the
% iteration through time, are only constructed once. Others are modified
% each time this function is called.

% The formulation and original code (for SISO systems) is by Elling W.
% Jacobsen from KTH university, Sweden. The formulation can be found here:
% http://www.kth.se/polopoly_fs/1.120161!/Menu/general/column-content/attachment/f13_ewj.pdf
% The code is modified and generalized for MIMO systems by Pooya Rezaei, 
% University of Vermont, USA. 

function [u,MPC_case] = MPC_calculation_delta_u(x0,yref_Np,MPC_case)

% MPC Set-up
nx = MPC_case.nx;
ny = MPC_case.ny;
nu = MPC_case.nu;
xref_vec = MPC_case.C_pinv*yref_Np;
xref_vec = [xref_vec; zeros(nu,MPC_case.Np)];
xref_vec = reshape(xref_vec,(nx+nu)*MPC_case.Np,1);

if ~isfield(MPC_case,'opt') % This condition is true only the first time this function is called
    % Define A_tilde, B_tilde, C_tilde for delta_u formulation
    MPC_case.A_tilde = [MPC_case.A, MPC_case.B; zeros(nu,nx), eye(nu)];
    MPC_case.B_tilde = [MPC_case.B; eye(nu)];
    MPC_case.C_tilde = [MPC_case.C, zeros(ny,nu)];
    
    % Compute mapping from present (initial) state, Ahat, and future inputs, Bhat,
    % to future states: x_future=Ahat*x0+Bhat*u_future
    Ahat = zeros(MPC_case.Np*(nx+nu),nx+nu);
    Bhat = zeros(MPC_case.Np*(nx+nu),MPC_case.Np*nu);
    for i = 1:MPC_case.Np
        idx_i = ((i-1)*(nx+nu)+1):i*(nx+nu);
        Ahat(idx_i,:) = MPC_case.A_tilde^i;
        for j = 1:MPC_case.Np
            idx_j = ((j-1)*nu+1):j*nu;
            if i >= j
                Bhat(idx_i,idx_j) = MPC_case.A_tilde^(i-j)*MPC_case.B_tilde;
            end
        end 
    end
    MPC_case.opt.Ahat = Ahat;
    MPC_case.opt.Bhat = Bhat;

    % Get system and MPC data
    MPC_case.Qx = MPC_case.C'*MPC_case.Qy*MPC_case.C;
    MPC_case.S = MPC_case.C'*MPC_case.Sy*MPC_case.C;
    % Stack up state and output weights in block-diagonal matrices to
    % replace summing in objective function by matrix multiplication in QP
    Qxhat = zeros(MPC_case.Np*(nx+nu),MPC_case.Np*(nx+nu));
    for i = 1:MPC_case.Np
        idx = (i-1)*(nx+nu)+1:(i-1)*(nx+nu)+nx;
        if i < MPC_case.Np
            Qxhat(idx,idx) = MPC_case.Qx;
        else
            Qxhat(idx,idx) = MPC_case.S;
        end
    end
%     Qxhat = sparse(Qxhat);
    MPC_case.opt.Qxhat = Qxhat;
    
    Quhat = zeros(nu*MPC_case.Np,nu*MPC_case.Np);
    for i = 1:MPC_case.Np
        idx = (i-1)*nu+1:i*nu;
        Quhat(idx,idx) = MPC_case.Qu;
    end        
%     Quhat = sparse(Quhat);
    
    % Constraints on states in matrix form: y_min<Hhat*x<y_max
    Hhat = zeros(MPC_case.Np*ny,MPC_case.Np*(nx+nu));
    for i = 1:MPC_case.Np
        idx_i = (i-1)*ny+1:i*ny;
        idx_j = (i-1)*(nx+nu)+1:i*(nx+nu);
        Hhat(idx_i,idx_j) = MPC_case.C_tilde;
    end        
%     Hhat = sparse(Hhat);
    MPC_case.opt.Hhat = Hhat;
    
    %% Sets up objective function and constraints for QP problem
    % Constraints
    Aineq_ymax = Hhat*Bhat;
    Aineq_ymin = -Aineq_ymax;
    Aineq_xmax = Bhat;
    Aineq_xmin = -Aineq_xmax;
    MPC_case.opt.Aineq = [Aineq_ymax; Aineq_ymin; Aineq_xmax; Aineq_xmin];

    % Lower & Upper bounds for variables
    MPC_case.opt.delta_u_max = repmat(MPC_case.delta_u_max,MPC_case.Np,1);
    MPC_case.opt.delta_u_min = repmat(MPC_case.delta_u_min,MPC_case.Np,1);

    % Obj func
    H = Bhat'*Qxhat*Bhat+Quhat;
    % H = (H+H')/2;
    round_digit = 9;
    MPC_case.opt.H = round(H*10^round_digit)/10^round_digit; % Added to avoid CPLEX error about H being not symmetric
end
Xdev = MPC_case.opt.Ahat*x0-xref_vec;
bineq_ymax = repmat(MPC_case.ymax,MPC_case.Np,1)-MPC_case.opt.Hhat*(Xdev+xref_vec);
bineq_ymin = repmat(-MPC_case.ymin,MPC_case.Np,1)+MPC_case.opt.Hhat*(Xdev+xref_vec);
bineq_xmax = repmat(MPC_case.xmax,MPC_case.Np,1)-(Xdev+xref_vec);
bineq_xmin = repmat(-MPC_case.xmin,MPC_case.Np,1)+(Xdev+xref_vec);

bineq = [bineq_ymax; bineq_ymin; bineq_xmax; bineq_xmin];
f = Xdev'*MPC_case.opt.Qxhat*MPC_case.opt.Bhat;

% Calculate future Np inputs:
[delta_u,~,exitflag] = quadprog(MPC_case.opt.H,f,MPC_case.opt.Aineq,bineq,[],[],MPC_case.opt.delta_u_min,MPC_case.opt.delta_u_max); 
x1 = MPC_case.A_tilde*x0+MPC_case.B_tilde*delta_u(1:nu);
u = x1(nx+1:nx+nu);
% cplexqp can also be used here with the same inputs and outputs if it is
% installed.
if exitflag<0
    disp('WARNING: infeasible QP problem.')
    pause(0.1)
end



