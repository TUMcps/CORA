function res = test_Krylov_intervalHomogeneousSolution(~)
% test_Krylov_intervalHomogeneousSolution - unit_test_function for checking
%    the Krylov method for the homegeneous solution of an initial interval
%
% Syntax:
%    res = test_Krylov_intervalHomogeneousSolution(~)
%
% Inputs:
%    no
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       23-August-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% system matrix
A = [...
  0.652541378199936  -0.292928516824949   0.549087928976002   0.274632260556064  -0.414286758065495  -0.141241992085791  -0.028985009015485   0.605937342611256  -0.195349385223348  -0.853531681405666; ...
   0.111507396391712   0.886771403763683  -0.926186402280864  -0.878814748052102  -0.240496674392058   0.630048246527933  -0.134643721501984  -0.098876392250540   0.534507822288518  -0.721762559941795; ...
  -0.184228225751299   0.728189979526127   0.844495081734710   0.095075956658037   0.414564407256731  -0.247618105941440   0.264265765751310  -0.006909214129143   -0.159428715921432   0.044005920605800; ...
   0.384842720632807  -0.427961712778828   0.159082783046048  -0.217631525890836   0.281950361270010  -0.782175454527015   0.295164577068571   0.228121531073434   0.360768593258595   0.343308641521092; ...
   0.466710516725629   0.647902449919577  -0.138568065795660  -0.048407199148869  -0.214935429197563  -0.371866410909703   0.794448727367410   0.886239821529420   0.230418230464495   0.487226934799377; ...
  -0.504441536402227  -0.683930638862332  -0.110954574785176   0.687389269618629   0.413797704097389   0.009111413582124   0.827318132599109   0.180993736522104   0.102327863083104  -0.176443263542827; ...
   0.015526934326904  -0.116581152633782   0.658696622614006  -0.174463584664272   0.205670047157241   0.622024548694704  -0.394835060811682  -0.452520224525244   0.388715632620385   0.060315005491066; ...
  -0.107197579271798  -0.080399191437818  -0.007555629564728  -0.262605107304710   0.522994075276099  -0.559668119751529  -0.817510097401123  -0.480638527682639   0.128110533548805  -0.041380923113904; ...
   0.268292332294290   0.714031800806371   0.395607487326902   0.400361847002023   0.525457042575656   0.216020969859562  -0.425914589519408  -0.031701286989668   -0.399781588176020  -0.805725638821021; ...
   0.216736942375567   0.516814156617457  -0.628881792412746  -0.110508717720818   0.514734156317458   0.112609184197350   0.197728711176581   0.022365142471784   0.486695871859436  -0.457708772484572];

% initial state
min = [...
   0.113338677026544; ...
  -0.364949053131162; ...
  -0.115799920838766; ...
  -1.130540023525459; ...
  -0.648439723643845; ...
  -0.674653320935817; ...
  -0.146071444296125; ...
  -1.366829171697049; ...
  -1.837093837409135; ...
  -0.267116720115053];

max = [...
   0.702468756588866; ...
   0.380322370568787; ...
   0.203903111033700; ...
  -0.400138900774135; ...
  -0.017619476598727; ...
   0.225921599407649; ...
   0.269245009861231; ...
  -1.080372480646494; ...
  -0.894093514119646; ...
   0.060086057402710];

x0 = interval(min,max);


%set options --------------------------------------------------------------
options.timeStep=0.1; %time step size for reachable set computation
options.tFinal = options.timeStep;
options.R0 = zonotope(x0);
options.x0 = center(x0);
options.U = zonotope(0*options.x0);
options.taylorTerms = 4;
options.krylovError = 1e-3;
options.reductionTechnique = 'girard';
options.zonotopeOrder = 1;
options.krylovStep = 1;
%--------------------------------------------------------------------------

%specify continuous dynamics-----------------------------------------------
linDyn = linearSys('KrylovTest',A,1); %initialize quadratic dynamics
linDyn_abs = linearSys('KrylovTest_abs',abs(A),1); %initialize quadratic dynamics
%--------------------------------------------------------------------------

% compute overapproximation------------------------------------------------
redOrder = ceil(0.5*length(options.x0));
b_c = center(x0);
b_delta = rad(x0);

% Arnoldi iteration
[V_c,H_c] = arnoldi(A, b_c, redOrder);
[V_delta,H_delta] = arnoldi(abs(A), b_delta, redOrder);

% solution without epsilon
expMatrix_delta = expm(H_delta*options.timeStep);
x_appr_delta = norm(b_delta)*V_delta*expMatrix_delta(:,1);

% compute epsilon values for A and |A|
epsilon = aux_computeEpsilon(linDyn, redOrder, options);
epsilon_abs = aux_computeEpsilon(linDyn_abs, redOrder, options);

% compute mu
mu = x_appr_delta + ones(length(b_c),1)*(norm(b_delta)*epsilon_abs*options.timeStep + norm(b_c)*epsilon*options.timeStep);

expMatrix_c = expm(H_c*options.timeStep);
x_appr = norm(b_c)*V_c*expMatrix_c(:,1) + interval(-mu, mu);
%--------------------------------------------------------------------------

% compute exact solution
x_exact = expm(A*options.timeStep)*zonotope(x0);

% Is exact solution in zonotope?
res = (interval(x_exact) <= x_appr);

end


% Auxiliary functions -----------------------------------------------------

function epsilon = aux_computeEpsilon(linSys, redOrder, options)
    %obtain time step, rho_A
    tau = options.timeStep;
    rho_A = normest(linSys.A);
    
    % auxiliary value for Wang
    [m,lambda,a] = aux_auxValues_Wang(linSys.A);
    q = aux_optimize_q_Wang(redOrder,tau,lambda,m);
    epsilon = aux_errorBound_Wang_normalized(rho_A,tau,q,a,lambda,redOrder);
end


function [m,lambda,a] = aux_auxValues_Wang(A)
% obtain a, b, and c
%     a = eigs((A+A')/2, 1, 'sm'); % lambda_min(A+A*/2)
%     b = eigs((A+A')/2, 1, 'lm'); % lambda_max(A+A*/2)
%     c = abs(eigs((A-A')/2, 1, 'lm')); % |lambda_max(A-A*/2)|
    a = eigs((A+A')/2, 1, 'sa'); % lambda_min(A+A*/2)
    b = eigs((A+A')/2, 1, 'la'); % lambda_max(A+A*/2)
    opts.tol = 1e-6;
    c = abs(eigs((A-A')/2, 1, 'li', opts)); % |lambda_max(A-A*/2)|

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

function q = aux_optimize_q_Wang(k,tau,lambda,m)
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

function err = aux_errorBound_Wang_normalized(rho_A,tau,q,a,lambda,KrylovOrder)
    
    % obtain remaining values
    Q = 11.08;
    exponent = -tau*(a-1/(2*lambda)*(1/q-q));
    err = 2*Q*tau*rho_A*q^(KrylovOrder-1)/(1-q)*exp(exponent); % Corr 5.3, Wang 2016
    
    %err = 2*Rinit_norm*Q*tau*rho_A*q^(KrylovOrder-1)/(1-q); % Corr 5.3, Wang 2016
end

% ------------------------------ END OF CODE ------------------------------
