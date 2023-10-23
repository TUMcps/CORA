function res = testMP_Krylov_linearTimeErrorBound(~)
% testMP_Krylov_linearTimeErrorBound - unit_test_function for checking
% whether the proposed error bound of the Krylov method is linear in time
% as derived in Appendix B of [1].
% This test requires the multiple precision toolbox.
%
% Syntax:
%    res = testMP_Krylov_linearTimeErrorBound(~)
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Reference:
%    [1] Althoff, M.: Reachability Analysis of Large Linear Systems With
%        Uncertain Inputs in the Krylov Subspace , IEEE Transactions on 
%        Automatic Control 65(2), 477–492, 2020.

% Authors:       Matthias Althoff
% Written:       17-November-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% enable access to private function "initReach_Krylov"

path = CORAROOT;
source = fullfile(path,'contDynamics','@linearSys','private','initReach_Krylov.m');
target = fullfile(path,'contDynamics','@linearSys','initReach_Krylov.m');
copyfile(source,target);
rmpath(genpath(path));
addpath(genpath(path));

% system with nu(A)<0 -----------------------------------------------------

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
options.R0 = zonotope([...
  -1.500147968675402  -0.919166483389253   0.740294075108247  -1.385060857214429  -0.128100444586379; ...
  -1.731247359232213  -1.515369480724521  -0.573387784727763  -0.989519038921135   0.279149822293821; ...
  -0.724963521890639   0.045574790470102  -0.194074587359824  -0.766919973323042  -0.999619885282059; ...
  -0.468234591594993   0.218855731178253  -1.620168750484433  -0.580708730517952  -1.859979909515807; ...
  -1.100108244156788  -0.944110286051283  -0.532136247490283  -1.264551181754348   0.669513515622473; ...
  -0.011762671798639  -0.414816979584528   0.156798725203194   0.120743459974948  -1.145060281809368; ...
  -0.995475582831655  -1.467319151693422  -0.118585500712028  -0.185247533976780  -0.893077985777429; ...
  -1.074866783228390  -0.310685806752072  -0.130252289367967  -0.032470808841968  -0.934033434879348; ...
  -0.980478681475524  -0.095671684159569  -0.942007809138542  -1.086728472520214   0.127262045423063; ...
   0.050513856730726  -0.109649516300879  -0.339706564980717  -1.451357086983043  -1.365878206307748]);


% set options
options.timeStep = 0.01; %time step size for reachable set computation
options.tFinal = 10*options.timeStep;
options.x0 = center(options.R0);
options.U = zonotope(0*options.x0);
options.uTrans = 0*options.x0;
options.taylorTerms = 4;
options.krylovError = 1e-6;
options.reductionTechnique = 'girard';
options.zonotopeOrder = 1;
options.krylovStep = 1;
%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %compute initial state factor
    options.factor(i)= options.timeStep^(i)/factorial(i);    
end

% specify continuous dynamics
linDyn = linearSys('KrylovTest',A,1); %initialize quadratic dynamics

% test errors of first system
res_1 = aux_testProcedure(linDyn, options);

% system with nu(A)>0 -----------------------------------------------------
N = 100;

% example matrix 1
i = 1:N;
diagVec = (i+1)./(N+1);
A = diag(diagVec);

% set options
options.timeStep = 0.1; %time step size for reachable set computation
options.tFinal = 10*options.timeStep;
options.x0 = exp(A)\ones(N,1);
options.R0 = zonotope(options.x0);
options.U = zonotope(0*options.x0);
options.uTrans = 0*options.x0;
options.taylorTerms = 4;
options.krylovError = 1e-6;
options.zonotopeOrder = 1;
options.krylovStep = 1;
%obtain factors for initial state and input solution
for i=1:(options.taylorTerms+1)
    %compute initial state factor
    options.factor(i)= options.timeStep^(i)/factorial(i);    
end

% specify continuous dynamics
linDyn = linearSys('KrylovTest',A,1); %initialize quadratic dynamics

% test errors of second system
res_2 = aux_testProcedure(linDyn, options);

% final result
res = res_1 && res_2;

% revoke access to private function "initReach_Krylov"
delete(target);
rmpath(genpath(path));
addpath(genpath(path));

end


% Auxiliary functions -----------------------------------------------------

function res = aux_testProcedure(linDyn, options)

% compute overapproximation
initReach_Krylov(linDyn, options.R0, options);

% obtain Krylov order
V = linDyn.krylov.state.V_c;
krylovOrder = length(V(1,:));

% check if order suffices accuracy for different times
% times
t = 0:options.timeStep:options.tFinal;

% check error for each time
res = [];
for i=1:length(t)
    % set error
    errorBound = options.krylovError*t(i);
    
    % compute Krylov error
    error = aux_KrylovError(linDyn,linDyn.A,options.x0,t(i),krylovOrder,options);
    
    % check whether error is within bounds
    res(end+1) = (error <= errorBound);
end

% All orders larger or equal?
res = all(res);

end


function error = aux_KrylovError(linDyn,A,v,t_f,krylovOrder,options)

    % set precision for variable precison toolbox
    precision = 100;

    % compute norm of v
    v_norm = norm(v);

    % initialize Krylov order and normalized error
    error = inf;
    nu_A = eigs((A+A')/2, 1, 'lm'); % lambda_max(A+A*/2) %<-- for Jia approx.
    nu_A = mp(nu_A, precision);
    dim = length(A);

    % Krylov order should not be larger than dimension
    if krylovOrder > dim
        krylovOrder = dim;
    end

    % times 
    delta = options.timeStep; % time step size
    if delta > t_f
        delta = t_f;
    end
    timeSteps = ceil(t_f/delta); % nr of time steps


    % perform Arnoldi iteration
    [~,H,Hlast,happyBreakdown] = arnoldi(A,v,krylovOrder+1);

    % init
    h_max = mp(0,precision);

    % compute error if happy breakdown did not occur
    if ~happyBreakdown
        for iStep=1:timeSteps

            %new time
            t = iStep*delta;

            %update values in Theorem 3.1 of Wang 2017
            % exponential matrix
            exp_tH = expm(H*t);
            if ~any(any(isnan(exp_tH))) 
                % extract entry
                h_t_Z = exp_tH(krylovOrder,1);
                % convert to interval matrix
                h_t_int = intervalMatrix(h_t_Z);
                % obtain absolute value
                h_t_tmp = abs(h_t_int.int);
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
        % determine phi
        phi = (exp(nu_A*t)-1)/nu_A;

        % error bound
        error = v_norm*Hlast*h_max*phi;

        if isnan(error) % if error is not a number (NaN)
            error = inf;
        end 
    end
end

% ------------------------------ END OF CODE ------------------------------
