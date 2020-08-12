function [V,H,krylovOrder] = subspace_Krylov_individual(A,nu_A,v,initKrylovOrder,options)
% subspace_Krylov_individual - computes the Krylov subspace for a single
% vector given the accuracy to be achieved; the a-posteriori approach in 
% equation 4.1 of 
%
% Jia, Z. & Lv, H. A posteriori error estimates of Krylov subspace 
% approximations to matrix functions Numerical Algorithms, 2015, 69, 1-28
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
% Written:      09-November-2018
% Last update:  13-November-2018          
% Last revision:---

%------------- BEGIN CODE --------------

% set precision for variable precison toolbox
precision = 34;

% compute norm of v
v_norm = norm(v);

% maximum allowed error
maxRelError = options.krylovError;

% initialize Krylov order and normalized error
krylovOrder = initKrylovOrder - 1;
errorBound_normalized = inf;
nu_A = mp(nu_A, precision);
dim = length(A);

% Krylov order should not be larger than dimension
if krylovOrder > dim
    krylovOrder = dim;
end

% times 
scalingFactor = 100;
t_f = options.tFinal; %final time
delta = scalingFactor*options.timeStep; % time step size
if delta > t_f
    delta = t_f;
end
timeSteps = ceil(t_f/delta); % nr of time steps

while (errorBound_normalized > maxRelError) && (krylovOrder <= dim)
    
    % increment Krylov order
    krylovOrder = krylovOrder + options.krylovStep;

    % perform Arnoldi iteration
    [V,H,Hlast,happyBreakdown] = arnoldi(A,v,krylovOrder+1);
    
%     % convert H to mp
%     H = mp(H,precision);
    
    % matrix zonotope of first time step
    Z_center = 0.5*H*delta;
    Z_delta{1} = Z_center;
    Htau = matZonotope(Z_center, Z_delta);
    
    % first exponential matrix
    maxOrder = 3;
    [exp_Htau_Z, exp_Htau_I] = expmOneParam(Htau,1,maxOrder);
    
    if any(any(isinf(supremum(exp_Htau_I.int))))
        %set error bound to infinity
        errorBound_normalized = inf;
    else

        % init
        h_max = mp(0,precision);

        % compute error if happy breakdown did not occur
        if ~happyBreakdown
            for iStep=1:timeSteps

                %new time
                t = (iStep-1)*delta;

                %update values in Theorem 3.1 of Wang 2017
                % exponential matrix
                exp_tH = expm(H*t);
                if ~any(any(isnan(exp_tH))) 
                    % exponential matrix for time interval
                    H_exp_Z = exp_tH*exp_Htau_Z;
                    H_exp_I = exp_tH*exp_Htau_I;
                    % extract entry
                    h_t_Z = H_exp_Z(krylovOrder,1);
                    h_t_I = H_exp_I(krylovOrder,1);
                    % convert to interval matrix
                    h_t_int = intervalMatrix(h_t_Z);
                    % obtain absolute value
                    h_t_tmp = abs(h_t_int.int+h_t_I);
                    h_t_abs = supremum(h_t_tmp);
                    % find maximum h_t
                    if h_t_abs > h_max
                        h_max = h_t_abs;
                    end
                else
                    % maximum value is set to infinity
                    h_max = inf;
                    break
                end
            end
            % determine phi_max
            if nu_A >= 0
                phi_max = (exp(nu_A*t_f)-1)/(nu_A*t_f);
            else
                phi_max = 1;
            end
            % error bound
            errorBound_normalized = v_norm*Hlast*h_max*phi_max;

            if isnan(errorBound_normalized) % if error is not a number (NaN)
                errorBound_normalized = inf;
            end 
        else
            % decrease Krylov order
            krylovOrder = length(V(1,:));
            break
        end
    end
end

%------------- END OF CODE --------------
