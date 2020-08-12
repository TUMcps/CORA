function [V,H,krylovOrder] = subspace_Krylov_individual_Wang(A,nu_A,v,initKrylovOrder,options)
% subspace_Krylov_individual - computes the Krylov subspace for a single
% vector given the accuracy to be achieved; the a-posteriori approach in 
% Theorem 3.1 of 
%
% Wang, H. & Ye, Q. Error Bounds for the Krylov Subspace Methods for 
% Computations of Matrix Exponentials SIAM Journal on Matrix Analysis and 
% Applications, 2017, 38, 155-187
%
% is used for tight error computation
%
% Syntax:  
%    [V,H] = subspace_Krylov_individual(obj,v,options)
%
% Inputs:
%    A - system matrix
%    nu_A - minimum eigenvalue of A
%    v - vector
%    initKrylovOrder - Krylov error that is first tested
%    options - reachability options
%
% Outputs:
%    V - orthonormal basis
%    H - Hessenberg matrix
%    KrylovOrder - dimension of the reduced system
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      06-November-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% set precision for variable precison toolbox
precision = 100;

% compute norm of v
v_norm = norm(v);

% maximum allowed error
maxRelError = options.krylovError;

% initialize Krylov order and normalized error
krylovOrder = initKrylovOrder;
errorBound_normalized = inf;
nu_A = mp(nu_A, precision);

% times 
t_f = options.tFinal; %final time
delta = options.timeStep; % time step size
timeSteps = ceil(t_f/delta); % nr of time steps


while errorBound_normalized > maxRelError

    % perform Arnoldi iteration
    [V,H,Hlast] = arnoldi(A,v,krylovOrder+1);

    % integrated h
    h_int = mp(0,precision);

    % init happy breakdown
    happyBreakdown = 0;

    for iStep=1:timeSteps

        %new time
        t = iStep*delta;

        %update values in Theorem 3.1 of Wang 2017
        H_exp = expm(-t*H);
        try
            h_t = H_exp(krylovOrder,1);
        catch
            %happy breakdown occurred
            happyBreakdown = 1;
            break
        end
        
        % add to h_int
        h_int = h_int + abs(h_t)*exp((t-t_f)*nu_A)*delta;
    end
    if ~happyBreakdown
        % error bound
        errorBound_normalized = v_norm*Hlast*h_int;

        if isnan(errorBound_normalized) % if error is not a number (NaN)
            errorBound_normalized = inf;
        end 
    else
        break
    end
    % increment Krylov order
    krylovOrder = krylovOrder + options.krylovStep;
end

%------------- END OF CODE --------------