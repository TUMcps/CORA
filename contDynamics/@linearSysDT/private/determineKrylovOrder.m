function [obj,options] = determineKrylovOrder(obj,options)
% determineKrylovOrder - determines the required order of the system in the
% Krylov subspace so that the maximum order is below the required threshold
%
% Syntax:  
%    options = determineKrylovOrder(obj,Rinit,options)
%
% Inputs:
%    obj - linearSys object
%    Rinit - initial set
%    options - options struct
%
% Outputs:
%    obj - linearSys object
%    options - options struct
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
% Written:      29-November-2016 
% Last update:  23-December-2016
%               28-October-2017
%               25-October-2018
% Last revision:---

%------------- BEGIN CODE --------------

% a priori bound
%[obj,options] = errorBound_Wang_aPriori(obj, options);

% a posteriori bound
[obj, options] = errorBound_Wang_aPosteriori(obj, options);
%[obj,options] = errorBound_Jia_aPosteriori(obj, options);

end

function err = errorBound_Saad(rho,Rinit_norm,KrylovOrder)
    err = 2*Rinit_norm*rho^KrylovOrder*exp(KrylovOrder)/factorial(KrylovOrder); % Thm 4.3, Saad 1992
    %err = Rinit_norm*rho^KrylovOrder/(2^(KrylovOrder-1)*factorial(KrylovOrder)); % eq (1), Stewart 1996
end

function [m,lambda,a] = auxValues_Wang(A)
    % obtain a, b, and c
%     a = eigs((A+A')/2, 1, 'sm'); % lambda_min(A+A*/2)
%     b = eigs((A+A')/2, 1, 'lm'); % lambda_max(A+A*/2)
%     c = abs(eigs((A-A')/2, 1, 'lm')); % |lambda_max(A-A*/2)|
%tic
    a = eigs((A+A')/2, 1, 'sa'); % lambda_min(A+A*/2)
    b = eigs((A+A')/2, 1, 'la'); % lambda_max(A+A*/2)
    opts.tol = 1e-6;
    c = abs(eigs((A-A')/2, 1, 'li', opts)); % |lambda_max(A-A*/2)|
%    toc

    % obtain alpha, beta from description after (4.1)
    alpha = (b-a)/2;
    beta = c;
    
    % equate lambda and lambda_1 to obtain m
    lambda_diff = inf;
    m_min = 0;
    m_max = 1;
    m = 0.5;
    while abs(lambda_diff) > eps
        % obtain m from (4.3)
        K = ellipticK(m); % elliptic integral of the first kind
        E = ellipticE(m); % elliptic integral of the second kind
        K_1 = ellipticK(1-m); % elliptic integral of the first kind
        E_1 = ellipticE(1-m); % elliptic integral of the second kind
        
        lambda = (E - (1-m)*K)/beta; % lambda is the ratio in (4.3)
        lambda_1 = (E_1 - m*K_1)/alpha; % lambda is the ratio in (4.3)
        lambda_diff = lambda - lambda_1;
        
        if lambda_diff < 0
            m_min = m;
            m = 0.5*(m + m_max);
        else
            m_max = m;
            m = 0.5*(m + m_min);
        end
    end
end

function q = optimize_q_Wang(k,tau,lambda,m)
    % init q
    q = 0.5;
    q_min = 0;
    q_max = 1;
    res = inf;
    
    % compute C
    C = tau/(2*lambda);
    
    while abs(res) > 1e-10
        res = (k-1)*q + (2-k)*q^2 - C*(1-q)*sqrt((1-q^2)^2 + 4*m*q^2);
        if res < 0
            q_min = q;
            q = 0.5*(q + q_max);
        else
            q_max = q;
            q = 0.5*(q + q_min);
        end
    end
end

function [err, WangExponent] = errorBound_Wang_normalized(rho_A,tau,q,a,m,lambda,KrylovOrder)
    % updated 05.12.2017 after paper of Wang appeared (better than arXiv version)
    % obtain remaining values
        
    zeta = 4*q^(KrylovOrder-1-tau*sqrt(m)/lambda)/(1-q)*min(tau*rho_A,q);
    L = 1/(2*lambda)*(1/q + q -2);
    sigma = exp(-tau*(a-L));
    err = zeta*max(sigma,1); %no division by tau since err is later multiplied by tau
    WangExponent = a-L; % the Wang exponent is used to determine whether the exponential matrix has to be recomputed
    
%     %STABLE MATRIX only:
%     q0 = 1/(a*lambda+1+sqrt((a*lambda+1)^2-1));
%     err = 4*q^(KrylovOrder-1-tau*sqrt(m)/lambda)/(1-q)*min(tau*rho_A,q0);
end

% function err = errorBound_Stewart(rho,Rinit_norm,lambda_min,KrylovOrder)
%     b = 2/(1+sqrt(5));
%     d = 1/(2+sqrt(5))*exp(b);
%     err = 4*Rinit_norm*exp(-lambda_min)*(exp(-b*KrylovOrder^2/rho)*(1+sqrt(rho*pi/(4*b))+d^rho/(1-d))); % eq (23), Stewart 1996; rho = \lambda_max - \lambda_min = \|A\|_2?
% end


function [obj,options] = errorBound_Wang_aPriori(obj, options)

    % maximum allowed error
    maxRelError = options.KrylovError;
    
    % auxiliary value for Wang
    [m,lambda,a] = auxValues_Wang(obj.A);
    rho_A = normest(obj.A);
    tau = options.timeStep;

    % matrix + time step
    A = obj.A*options.timeStep;

    % compute system matrix norm
    rho = normest(A);

    % % smallest eigenvalue (real part or absolute?)
    % lambda_min = abs(eigs(A, 1, 'sm'));
    % lambda_max = abs(eigs(A, 1, 'lm'));
    % rho_Stewart = lambda_max - lambda_min;

    % init 
    KrylovOrder = 0;
    errorBound_normalized = inf;

    % increase reduced order step by step
    while errorBound_normalized > maxRelError
        % increment Krylov order
        KrylovOrder = KrylovOrder +1;

        q = optimize_q_Wang(KrylovOrder,tau,lambda,m);

        %q = 0.1;

        % relative error bound
        [errorBound_normalized, WangExponent] = errorBound_Wang_normalized(rho_A,tau,q,a,m,lambda,KrylovOrder);
    end


    % save Krylov order in options struct
    options.redDim = KrylovOrder;

    % save results to object structure
    obj.krylov.Atnorm = rho;
    obj.krylov.errorBound_normalized = errorBound_normalized;
    obj.krylov.WangExponent = WangExponent;
end

function [obj,options] = errorBound_Wang_aPosteriori(obj, options)
    % 01.11.2018
    
    % change obj.A to fit paper
    A = -obj.A;
    
    % set precision for variable precison toolbox
    precision = 100;
    
    % maximum allowed error
    maxRelError = options.KrylovError;
    
    % initialize Krylov order and normalized error
    KrylovOrder = 1;
    errorBound_normalized = inf;
    
    % dimension
    dim = length(obj.A);
    
    % times 
    t_f = options.tFinal; %final time
    delta = options.timeStep; % time step size
    timeSteps = ceil(t_f/delta); % nr of time steps
    
    % initialize initial vector
    vInit = [1;zeros(dim-1,1)];
    
    % unit vector of first coordinate
    e_1 = [1; zeros(dim-1,1)];
    
    % minimum eigenvalue
    nu_A = eigs((A+A')/2, 1, 'sa'); % lambda_min(A+A*/2)
    nu_A = mp(nu_A, precision);
    
    while errorBound_normalized > maxRelError
        % increment Krylov order
        %KrylovOrder = KrylovOrder + 1;
        KrylovOrder = KrylovOrder + 10;

        % perform Arnoldi iteration
        [~,H,Hlast] = arnoldi(A,vInit,KrylovOrder+1);

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
                h_t = H_exp(KrylovOrder,1);
            catch
                %happy breakdown occurred
                happyBreakdown = 1;
                break
            end

%             % exponent correction
%             exp_corr = log(abs(h_t));
% 
%             % add to h_int
%             h_int = h_int + exp((t-t_f)*nu_A+exp_corr)*delta;

            % add to h_int
            h_int = h_int + abs(h_t)*exp((t-t_f)*nu_A)*delta;
        end
        if ~happyBreakdown
            % error bound
            errorBound_normalized = Hlast*h_int;

            if isnan(errorBound_normalized) % if error is not a number (NaN)
                errorBound_normalized = inf;
            end 
            KrylovOrder
            H(KrylovOrder+1,KrylovOrder)
            errorBound_normalized
        else
            break
        end
    end
    
    % save Krylov order in options struct
    if happyBreakdown
        % error bound
        options.redDim = length(H_exp(:,1));
        obj.krylov.errorBound_normalized_tf = 0;
    else     
        options.redDim = KrylovOrder;
        obj.krylov.errorBound_normalized_tf = errorBound_normalized;
    end
    
end

function [obj,options] = errorBound_Jia_aPosteriori(obj, options)
    % 01.11.2018
    
    % change obj.A to fit paper
    A = -obj.A;
    
    % set precision for variable precison toolbox
    precision = 100;
    
    % maximum allowed error
    maxRelError = options.KrylovError;
    
    % initialize Krylov order and normalized error
    KrylovOrder = 1;
    errorBound_normalized = inf;
    
    % dimension
    dim = length(obj.A);
    
    % times 
    t_f = options.tFinal; %final time
    delta = options.timeStep; % time step size
    timeSteps = ceil(t_f/delta); % nr of time steps
    
    % initialize initial vector
    vInit = [1;zeros(dim-1,1)];
    
    % unit vector of first coordinate
    e_1 = [1; zeros(dim-1,1)];
    
    % minimum eigenvalue
    nu_A = eigs((A+A')/2, 1, 'lm'); % lambda_max(A+A*/2)
    nu_A = mp(nu_A, precision);
    
    while errorBound_normalized > maxRelError
        % increment Krylov order
        %KrylovOrder = KrylovOrder + 1;
        KrylovOrder = KrylovOrder + 10;

        % perform Arnoldi iteration
        [~,H,Hlast] = arnoldi(A,vInit,KrylovOrder+1);
    
        % init
        h_max = mp(0,precision);

        for iStep=1:timeSteps
            
            %new time
            t = iStep*delta;

            %update values in Theorem 3.1 of Wang 2017
            H_exp = expm(-t*H);
            try
                h_t = H_exp(KrylovOrder,1);
            catch
                %happy breakdown occurred
                happyBreakdown = 1;
                break
            end
            
            % find maximum h_t
            if h_t > h_max
                h_max = h_t;
            end
        end
        % exponent correction
        %exp_corr = log(abs(h_t));
        
        % error bound
        errorBound_normalized = Hlast*h_max*(exp(-t_f*nu_A)-1)/(-nu_A);
        if isnan(errorBound_normalized) % if error is not a number (NaN)
            errorBound_normalized = inf;
        end 
        KrylovOrder
        H(KrylovOrder+1,KrylovOrder)
        errorBound_normalized
    end
    % save Krylov order in options struct
    options.redDim = KrylovOrder;
    obj.krylov.errorBound_normalized_tf = errorBound_normalized;
end

%------------- END OF CODE --------------