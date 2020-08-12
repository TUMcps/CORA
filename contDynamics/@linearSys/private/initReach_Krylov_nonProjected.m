function [Rfirst,options] = initReach_Krylov(obj,Rinit,options)
% reach - computes the reachable continuous set for the first time step
% using Krylov subspace methods
%
% Syntax:  
%    [obj,Rfirst,options] = initReach_Krylov(obj,Rinit,options)
%
% Inputs:
%    obj - linearSys object
%    Rinit - initial reachable set
%    options - options for the computation of the reachable set
%
% Outputs:
%    obj - linearSys object
%    Rfirst - first reachable set 
%    options - options for the computation of the reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      15-November-2016
% Last update:  29-November-2016
%               19-December-2016
%               03-March-2017
%               02-November-2018
% Last revision:---

%------------- BEGIN CODE --------------

% obtain reduced dimension from demanded maximum error
[obj,options] = determineKrylovOrder(obj,options);

% INITIAL SET SOLUTION------------------------------------------------------
% obtain center and generators of initial set
c = sparse(center(Rinit));
G = sparse(generators(Rinit));

% check if center is zero
c_norm = norm(c);
if c_norm == 0
    c_new = c;
    tie_error = 0*c;
else
    %Arnoldi
    [V_c,H_c] = arnoldi(obj,c,options.redDim);
    
    % generate reduced order system for the center
    A_red = H_c;
    B_red = V_c'*obj.B;
    linRedSys_c = linearSys('linearReducedDynamics',A_red,B_red); %initialize linear reduced dynamics
    
    % compute exponential matrix
    linRedSys_c = exponential(linRedSys_c,options);
    
    % compute time interval error (tie)
    linRedSys_c = tie(linRedSys_c,options);
    F_tmp = linRedSys_c.taylor.F;
    tie_error = V_c*F_tmp*V_c.'*c;
    
    %Compute new center
    eAt = expm(A_red*options.timeStep);
    c_new = c_norm*V_c*eAt(:,1);
    
    %store error
    %taylor_error_c = linRed_c.taylor.error; %required for input solution
end


% preallocation
dim = length(c);
nrOfGens = length(G(1,:));
g_norm = zeros(nrOfGens,1);
V_g = cell(nrOfGens);
H_g = cell(nrOfGens);
G_new = zeros(dim,nrOfGens);

% obtain generators using the Arnoldi iteration
for iGen = 1:nrOfGens
    g_norm(iGen) = norm(G(:,iGen));
    if g_norm(iGen) == 0
        G_new(:,iGen) = G(:,iGen);
    else
        %Arnoldi
        [V_g{iGen},H_g{iGen}] = arnoldi(obj,G(:,iGen),options.redDim);
        
        % generate reduced order system for the generators
        A_red_g = H_g{iGen};
        B_red_g = V_g{iGen}'*obj.B;
        linRedSys_g = linearSys('linearReducedDynamics',A_red_g,B_red_g); %initialize linear reduced dynamics
        
        % compute exponential matrix
        linRedSys_g = exponential(linRedSys_g,options);
        
        % compute time interval error (tie)
        linRedSys_g = tie(linRedSys_g,options);
        F_tmp = linRedSys_g.taylor.F;
        %tie_error 
        new_tie_error_mid = V_g{iGen}*center(F_tmp)*V_g{iGen}.'*G(:,iGen);
        new_tie_error_rad = abs(V_g{iGen}*rad(F_tmp)*V_g{iGen}.'*G(:,iGen));
        tie_error = tie_error + interval(new_tie_error_mid - new_tie_error_rad, new_tie_error_mid + new_tie_error_rad);
        
        %Compute new generator
        eAt = expm(A_red_g*options.timeStep);
        G_new(:,iGen) = g_norm(iGen)*V_g{iGen}*eAt(:,1);
        
        %store error
        %taylor_error_g{iGen} = linRed_g{iGen}.taylor.error; %required for input solution
    end
    %iGen
end

% Krylov error computation
if options.krylovError > 2*eps
    error_normalized = obj.krylov.errorBound_normalized;
    error = norm(Rinit)*error_normalized;
    Krylov_interval = interval(-1,1)*ones(dimension(obj),1)*error;
    R_Krylov = zonotope(Krylov_interval);
    
    % initial-state-solution zonotope
    R_initSol = zonotope([c_new,G_new,error*eye(dimension(obj))]);
else
    R_Krylov = 0;
    
    % initial-state-solution zonotope
    R_initSol = zonotope([c_new,G_new]);
end

% tie zonotope
R_tie = zonotope(tie_error);

%--------------------------------------------------------------------------

% INPUT SOLUTION-----------------------------------------------------------
obj = inputSolution_Krylov(obj, options);
%--------------------------------------------------------------------------

% change the time step
linRed.taylor.timeStep = options.timeStep;

% load data
Rtrans = obj.taylor.Rtrans;
RV = obj.taylor.RV;
inputCorr = obj.taylor.inputCorr;

% first time step homogeneous solution
Rhom_tp = R_initSol + R_Krylov + Rtrans; %Rtrans not yet changed
Rhom = enclose(Rinit,Rhom_tp) + R_tie + RV + inputCorr;

%reduce zonotope
Rhom=reduce(Rhom,options.reductionTechnique,options.zonotopeOrder);
Rhom_tp=reduce(Rhom_tp,options.reductionTechnique,options.zonotopeOrder);
RV=reduce(RV,options.reductionTechnique,options.zonotopeOrder);

%save homogeneous and particulate solution
options.Rhom=Rhom;
options.Rhom_tp=Rhom_tp;
options.Raux=RV;
options.Rpar=interval(RV);
options.Rtrans=obj.taylor.Rtrans;
options.R_tie = R_tie;
options.inputCorr = inputCorr;
options.Rtrans = Rtrans;
options.Rhom_0 = Rhom;
options.RV_0 = RV;


%total solution
if isa(Rinit,'mptPolytope')
    %convert zonotopes to polytopes
    Radd = mptPolytope(RV);
    Rtotal = Rhom+Radd;
    Rtotal_tp = Rhom_tp+Radd;
else
    %original computation
    Rtotal = Rhom+RV;
    Rtotal_tp = Rhom_tp+RV;
end

%write results to reachable set struct Rfirst
Rfirst.tp=Rtotal_tp;
Rfirst.ti=Rtotal;


%precompute solutions for state solution
[V_c,H_c,V_g,H_g] = subspace_Krylov(obj,Rinit,options);

% save precomputed solutions
obj.krylov.state.V_c = V_c;
obj.krylov.state.H_c = H_c;
obj.krylov.state.V_g = V_g;
obj.krylov.state.H_g = H_g;

%precompute solutions for input solution
[V_c,H_c,V_g,H_g] = subspace_Krylov(obj,RV,options);

% save precomputed solutions
obj.krylov.input.V_c = V_c;
obj.krylov.input.H_c = H_c;
obj.krylov.input.V_g = V_g;
obj.krylov.input.H_g = H_g;



%------------- END OF CODE --------------