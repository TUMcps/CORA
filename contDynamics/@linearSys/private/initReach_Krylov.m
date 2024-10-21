function [Rfirst,options] = initReach_Krylov(linsys,Rinit,options)
% initReach_Krylov - computes the reachable continuous set for the first
%    time step using Krylov subspace methods
%
% Syntax:
%    [Rfirst,options] = initReach_Krylov(linsys,Rinit,options)
%
% Inputs:
%    linsys - linearSys object
%    Rinit - initial reachable set
%    options - options for the computation of the reachable set
%
% Outputs:
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

% Authors:       Matthias Althoff
% Written:       15-November-2016
% Last update:   29-November-2016
%                19-December-2016
%                03-March-2017
%                02-November-2018
%                24-July-2020 (box inputs removed)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% retrieve C matrix for projection
C = linsys.C;

% precompute solutions for state solution
[V_c,H_c,V_g,H_g] = subspace_Krylov(linsys,Rinit,options);

% save precomputed solutions and project them
% center
if ~isempty(V_c)
    linsys.krylov.state.V_c = V_c;
    linsys.krylov.state.V_c_proj = C*V_c;
else
    linsys.krylov.state.V_c_proj = [];
end
%generators
if ~isempty(V_g)
    for i=1:length(V_g)
        linsys.krylov.state.V_g{i} = V_g{i};
        linsys.krylov.state.V_g_proj{i} = C*V_g{i};
    end
else
    linsys.krylov.state.V_g_proj = [];
end
linsys.krylov.state.H_c = H_c;
linsys.krylov.state.H_g = H_g;


% INITIAL SET SOLUTION------------------------------------------------------
% obtain center and generators of initial set
c = sparse(Rinit.c);
G = sparse(Rinit.G);

% check if center is zero
c_norm = norm(c);
if c_norm == 0
    c_new = c;
    tie_error = 0*c;
else
    % generate reduced order system for the center
    A_red = linsys.krylov.state.H_c;
    B_red = linsys.krylov.state.V_c'*linsys.B;
    % initialize linear reduced dynamics
    linsys_red_c = linearSys('linearReducedDynamics',A_red,B_red);
    
    % compute exponential matrix and correction matrix for the state
    [~,F] = taylorMatrices(linsys_red_c,options.timeStep,options.taylorTerms);
    
    % compute time interval error (tie)
    Fmid = center(F);
    Frad = rad(F);
    norm_c = norm(c);
    new_tie_error_mid = norm_c*V_c*Fmid(:,1);
    new_tie_error_rad = norm_c*abs(V_c*Frad(:,1));
    tie_error = interval(new_tie_error_mid - new_tie_error_rad, ...
        new_tie_error_mid + new_tie_error_rad);

    
    %Compute new center
    eAt = expm(A_red*options.timeStep);
    c_new = c_norm*V_c*eAt(:,1);
    
    %store error
    %taylor_error_c = linRed_c.taylor.error; %required for input solution
end


% preallocation
n = length(c);
nrOfGens = size(G,2);
g_norm = zeros(nrOfGens,1);
% V_g = obj.krylov.state.V_g;
% H_g = obj.krylov.state.H_g;
G_new = zeros(n,nrOfGens);

% obtain generators using the Arnoldi iteration
for iGen = 1:nrOfGens
    g_norm(iGen) = norm(G(:,iGen));
    if g_norm(iGen) == 0
        G_new(:,iGen) = G(:,iGen);
    else

        % generate reduced order system for the generators
        A_red_g = H_g{iGen};
        B_red_g = V_g{iGen}'*linsys.B;
        %initialize linear reduced dynamics
        linsys_red_g = linearSys('linearReducedDynamics',A_red_g,B_red_g);
        
        % compute exponential matrix
        [~,F] = taylorMatrices(linsys_red_g,options.timeStep,options.taylorTerms);
        
        % compute time interval error (tie)
        Fmid = center(F);
        Frad = rad(F);
        norm_g = norm(G(:,iGen));
        new_tie_error_mid = norm_g*V_g{iGen}*Fmid(:,1);
        new_tie_error_rad = norm_g*abs(V_g{iGen}*Frad(:,1));
        tie_error = tie_error + ...
            interval(new_tie_error_mid - new_tie_error_rad, ...
                     new_tie_error_mid + new_tie_error_rad);
        
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
%     error_normalized = obj.krylov.errorBound_normalized;
%     error = norm(Rinit)*error_normalized;
    err = options.krylovError;
    Krylov_interval = interval(-1,1)*ones(linsys.nrOfStates,1)*err;
    R_Krylov = zonotope(Krylov_interval);
    
    % initial-state-solution zonotope
    R_initSol = zonotope(c_new,G_new,err*eye(linsys.nrOfStates));
else
    R_Krylov = 0;
    
    % initial-state-solution zonotope
    R_initSol = zonotope(c_new,G_new);
end

% tie zonotope
R_tie = zonotope(tie_error);

%--------------------------------------------------------------------------

% INPUT SOLUTION-----------------------------------------------------------
linsys = inputSolution_Krylov(linsys, options);
%--------------------------------------------------------------------------

% load data
Rtrans = linsys.taylor.Rtrans;
RV = linsys.taylor.RV;
RV_0 = RV + Rtrans;
C_trans = C;
inputCorr = linsys.taylor.inputCorr;

% first time step homogeneous solution
Rhom_tp = R_initSol + R_Krylov + Rtrans; %Rtrans not yet changed
Rhom = enclose(Rinit,Rhom_tp) + R_tie + RV + inputCorr;

%reduce zonotope
Rhom=reduce(Rhom,options.reductionTechnique,options.zonotopeOrder);
Rhom_tp=reduce(Rhom_tp,options.reductionTechnique,options.zonotopeOrder);
RV=reduce(RV,options.reductionTechnique,options.zonotopeOrder);

% precompute solutions for input solution 
[V_c,H_c,V_g,H_g] = subspace_Krylov(linsys,RV_0,options);

% save precomputed solutions and project them
% center 
if ~isempty(V_c)
    linsys.krylov.input.V_c_proj = C_trans*V_c;
else
    linsys.krylov.input.V_c_proj = [];
end
% generators
linsys.krylov.input.V_g_proj = []; % initialize so that field 'V_g_proj' always exists
for i=1:length(V_g)
    if ~isempty(V_g{i})
        linsys.krylov.input.V_g_proj{i} = C_trans*V_g{i};
    else
        linsys.krylov.input.V_g_proj{i} = [];
    end
end

linsys.krylov.input.H_c = H_c;
linsys.krylov.input.H_g = H_g;

%save homogeneous and particulate solution
options.Rhom_proj = C*Rhom;
options.Rhom_tp_proj = C*Rhom_tp;
options.Raux_proj = C*RV;
options.Rpar_proj = interval(C*RV);
options.Rtrans_proj = C*linsys.taylor.Rtrans;
options.R_tie_proj = C*R_tie;
options.inputCorr_proj = C*inputCorr;
options.Rhom_0 = Rhom;
options.RV_0 = RV_0;

%total solution
if isa(Rinit,'polytope')
    %convert zonotopes to polytopes
    Radd = polytope(RV);
    Rtotal = Rhom+Radd;
    Rtotal_tp = Rhom_tp+Radd;
else
    %original computation
    Rtotal = Rhom+RV;
    Rtotal_tp = Rhom_tp+RV;
end

%write results to reachable set struct Rfirst
Rfirst.tp = C*Rtotal_tp;
Rfirst.ti = C*Rtotal;

% ------------------------------ END OF CODE ------------------------------
