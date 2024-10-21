function sys = platoonN(h,n_v,n_k)
% platoonN - system dynamics for the discrete-time version of the platoon 
%   benchmark with n_v vehicles (see Sec. VII in [1])
%
% Syntax:  
%    sys = platoonN(x,u,h)
%
% Inputs:
%    h - time step size
%    N_v - number of follower vehicles
%
% Outputs:
%    sys - time-varying linearSysDT
% 
% References:
%    [1] M. Althoff et al. "Reachability analysis of nonlinear systems with 
%        uncertain parameters using conservative linearization", in Proc. 
%        of the 62nd IEEE Conference on Decision and Control, pp.
%        4042-4048, 2008.

% Author:        Laura Luetzow, Matthias Althoff
% Written:       31-January-2024
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

A = cell(n_k,1);
B = cell(n_k,1);
C = cell(n_k,1);
D = cell(n_k,1);

for k = 1:n_k
    % time-varying parameter
    xi = 2 + 0.3*sin(0.2*k*h);

    % differential equations
    B_cont = [];
    B_column = [xi; 0; 0];
    A_0 = [-xi 0 0; 0 0 1; -1 0 0];
    A_i = [zeros(3,3) A_0];
    A_i(3,1) = 1;
    for i = 0:n_v-1
        if i==0
            A_cont = [A_0 zeros(3, 3*(n_v-1))];
        else
            A_cont = [A_cont; zeros(3,(i-1)*3) A_i zeros(3,(n_v-i-1)*3)];
        end
        B_cont = [B_cont; zeros(3,i) B_column zeros(3,n_v-i-1)];
    end
    A{k} = expm(A_cont*h);
    B{k} = integral(@(t)expm(A_cont*t),0,h,'ArrayValued',true)*B_cont;
    C{k} = eye(length(A{k}));
    D{k} = [];
end

sys = linearSysDT(sprintf('platoon%d',n_v),A,B,[],C,D,h);
end

%------------- END OF CODE --------------
