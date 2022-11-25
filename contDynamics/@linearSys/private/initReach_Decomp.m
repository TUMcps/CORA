function [Yhat0,options] = initReach_Decomp(obj,Rinit,options)
% initReach_Decomp - computes the reachable continuous output set
%  for the first time step (time interval and time point)
%
% Syntax:  
%    [Yhat0,options] = initReach_Decomp(obj,Rinit,options)
%
% Inputs:
%    obj     - linearSys object
%    Rinit   - initial reachable set
%    options - options for the computation of the reachable set
%
% Outputs:
%    Rfirst  - first reachable set 
%    options - options for the computation of the reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      28-May-2019 
% Last update:  09-Aug-2019
% Last revision:---

%------------- BEGIN CODE --------------


% decide sparsity ---------------------------------------------------------
sparseMatrices = false;
if issparse(obj.A)
    sparseMatrices = true;
elseif obj.dim >= 10 % dim at least 10 for sparsity to be useful 
    % A non-sparse data-structure, but conversion could be useful
    logicalA = obj.A ~= 0;
    density = sum(sum(logicalA)) / obj.dim^2;
    if density < 0.5 % ... adapt this factor if necessary
        sparseMatrices = true;
    end
end
% -------------------------------------------------------------------------


% obtain factors for initial state and input solution ---------------------
% not used in exponential, but in inputTie
options.factor = zeros(options.taylorTerms+1,1);
for i=1:(options.taylorTerms+1)
    % compute initial state factor
    options.factor(i) = options.timeStep^(i)/factorial(i);    
end
% -------------------------------------------------------------------------


% preliminary computations ------------------------------------------------
% compute exponential matrix
obj = exponential(obj,options);
% compute time interval error (tie)
obj = tie(obj,options);
% compute reachable set due to input
obj = inputSolution(obj,options);
% save the time step
obj.taylor.timeStep = options.timeStep;

% compute reachable set of first time interval
if sparseMatrices
    eAt = sparse(expm(obj.A*options.timeStep));
else
    eAt = expm(obj.A*options.timeStep);
end
% save data to object structure
obj.taylor.eAt = eAt;
% -------------------------------------------------------------------------


% read sets from preliminary computations ---------------------------------
F = obj.taylor.F;
RV = obj.taylor.RV;
inputCorr = obj.taylor.inputCorr;
if iscell(obj.taylor.Rtrans)
    Rtrans = obj.taylor.Rtrans{1};
else
    Rtrans = obj.taylor.Rtrans;
end
% -------------------------------------------------------------------------


% first time step homogeneous solution ------------------------------------
Rhom_tp = eAt*Rinit + Rtrans;
if isa(Rinit,'quadZonotope') 
    Rhom = enclose(Rinit,Rhom_tp) + F*zonotope(Rinit) + inputCorr;
elseif isa(Rinit,'zonoBundle') 
    Rhom = enclose(Rinit,Rhom_tp) + F*Rinit.Z{1} + inputCorr;
else
    Rhom = enclose(Rinit,Rhom_tp) + F*Rinit + inputCorr;
end
% -------------------------------------------------------------------------


% reduce zonotopes --------------------------------------------------------
Rhom = reduce(Rhom,options.reductionTechnique,options.zonotopeOrder);
Rhom_tp = reduce(Rhom_tp,options.reductionTechnique,options.zonotopeOrder);
RV = reduce(RV,options.reductionTechnique,options.zonotopeOrder);
% -------------------------------------------------------------------------


%save homogeneous and particulate solution --------------------------------
options.Rhom = Rhom;
options.Rhom_tp = Rhom_tp;
options.Raux = RV;
options.Rpar = RV;
options.Rtrans = obj.taylor.Rtrans;
% -------------------------------------------------------------------------


%total solution -----------------------------------------------------------
if isa(Rinit,'mptPolytope')
    %convert zonotopes to polytopes
    Radd = mptPolytope(RV);
    Rtotal = Rhom + Radd;
    Rtotal_tp = Rhom_tp + Radd;
else
    %original computation
    Rtotal = Rhom+RV;
    Rtotal_tp = Rhom_tp+RV;
end
% -------------------------------------------------------------------------


% decompose solution to blocks --------------------------------------------
Xhat0.ti = cell(options.blocks,1); % -> Yhat0.ti
Xhat0.tp = cell(options.blocks,1); % -> Yhat0.tp
Yhat0.ti = cell(options.blocks,1); % -> Rout
Yhat0.tp = cell(options.blocks,1); % -> Rout_tp
C        = cell(options.blocks,1);
Cno0     = false(options.blocks,1);
% -------------------------------------------------------------------------


% containers for each iteration -------------------------------------------
Rhom0    = cell(options.blocks,1); % -> options.Rhom in blocks
Rhom_tp0 = cell(options.blocks,1); % -> options.Rhom_tp in blocks
Rinhom   = cell(options.blocks,1); % -> P^R in blocks
Yinhom   = cell(options.blocks,1); % -> C * P^R in blocks
Rtrans   = cell(options.blocks,1); % -> options.Rtrans in blocks
Raux     = cell(options.blocks,1); % -> options.Raux in blocks
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
for b=1:options.blocks
    s_i = options.partition(b,1); % start of row-block
    f_i = options.partition(b,2); % end of row-block
    
    C{b} = obj.C(:,s_i:f_i);
    Cno0(b) = nnz(C{b}) ~= 0; % if b-th block of C is not all zeros
    
    Rfirst.ti{b,1} = project(Rtotal,s_i:f_i);
    Rfirst.tp{b,1} = project(Rtotal_tp,s_i:f_i);
    
    % first step without integration of C into X for calculation of Y
    if sparseMatrices
        Xhat0.ti{b} = zonotope([center(Rfirst.ti{b}),sparse(generators(Rfirst.ti{b}))]);
        Xhat0.tp{b} = zonotope([center(Rfirst.tp{b}),sparse(generators(Rfirst.tp{b}))]);
        % decompose homogeneous solution
        Rhom0{b} = zonotope([center(project(options.Rhom,s_i:f_i)),...
            sparse(generators(project(options.Rhom,s_i:f_i)))]);
        Rhom_tp0{b} = zonotope([center(project(options.Rhom_tp,s_i:f_i)),...
            sparse(generators(project(options.Rhom_tp,s_i:f_i)))]);
    else
        Xhat0.ti{b} = zonotope([center(Rfirst.ti{b}),generators(Rfirst.ti{b})]);
        Xhat0.tp{b} = zonotope([center(Rfirst.tp{b}),generators(Rfirst.tp{b})]);
        % decompose homogeneous solution
        Rhom0{b} = zonotope([center(project(options.Rhom,s_i:f_i)),...
            generators(project(options.Rhom,s_i:f_i))]);
        Rhom_tp0{b} = zonotope([center(project(options.Rhom_tp,s_i:f_i)),...
            generators(project(options.Rhom_tp,s_i:f_i))]);
    end  
    
end
% -------------------------------------------------------------------------

% simplification if no inhomogenuity (skip computation)
options.isInhom = nnz(options.uTrans) ~= 0 || nnz(generators(options.U)) ~= 0;

% -------------------------------------------------------------------------
for bi=1:options.blocks
    % i - row start and end of block
    s_i = options.partition(bi,1);
    f_i = options.partition(bi,2); 

	% take row block of Raux and Rtrans -----------------------------------
    if sparseMatrices
        Raux{bi} = zonotope([center(project(options.Raux,s_i:f_i)),...
            sparse(generators(project(options.Raux,s_i:f_i)))]);
        Rtrans{bi} = zonotope([center(project(options.Rtrans,s_i:f_i)),...
            sparse(generators(project(options.Rtrans,s_i:f_i)))]);
    else
        Raux{bi} = zonotope([center(project(options.Raux,s_i:f_i)),...
            generators(project(options.Raux,s_i:f_i))]);
        Rtrans{bi} = zonotope([center(project(options.Rtrans,s_i:f_i)),...
            generators(project(options.Rtrans,s_i:f_i))]);
    end
    % ---------------------------------------------------------------------

    
    % calculate output blocks (only if corresponding C block non-zero -----
    if Cno0(bi)
        Yhat0.ti{bi} = full(C{bi} * Xhat0.ti{bi});
        Yhat0.tp{bi} = full(C{bi} * Xhat0.tp{bi});
        
        % calculate inhomogeneous solutions for state + output ------------
        if options.isInhom

            % no propagation of options.Rtrans in init step, V ~ options.Raux
            if sparseMatrices
                Rinhomtemp = zonotope([center(project(options.Rtrans,s_i:f_i)),...
                    sparse(generators(project(options.Rpar,s_i:f_i) + project(options.Rtrans,s_i:f_i)))]);
                Yinhomtemp = zonotope([center(C{bi}*project(options.Rtrans,s_i:f_i)),...
                    sparse(C{bi}*generators(project(options.Rpar,s_i:f_i) + project(options.Rtrans,s_i:f_i)))]);
            else
                Rinhomtemp = zonotope([center(project(options.Rtrans,s_i:f_i)),...
                    generators(project(options.Rpar,s_i:f_i) + project(options.Rtrans,s_i:f_i))]);
                Yinhomtemp = zonotope([center(C{bi}*project(options.Rtrans,s_i:f_i)),...
                    C{bi}*generators(project(options.Rpar,s_i:f_i) + project(options.Rtrans,s_i:f_i))]);
            end

            % propagate Rinhomtemp / Yinhomtemp over all blocks
            for bj=1:options.blocks
                % j - col start and end of block
                s_j = options.partition(bj,1);
                f_j = options.partition(bj,2);
                eAt_block = eAt(s_i:f_i,s_j:f_j);

                if nnz(eAt_block)
                    % don't enter if block all-zero
                    % both Rtrans and Raux needed, Rinhom = V (incl. in 2nd zono)
                    if sparseMatrices
                        Rinhomtemp = Rinhomtemp + ...
                            zonotope([eAt_block * center(project(options.Raux,s_j:f_j)),...
                            sparse(eAt_block * generators(project(options.Raux,s_j:f_j)))]);
                        Yinhomtemp = Yinhomtemp + ...
                            zonotope([(C{bi} * eAt_block) * center(project(options.Raux,s_j:f_j)),...
                            sparse((C{bi} * eAt_block) * generators(project(options.Raux,s_j:f_j)))]);
                    else
                        Rinhomtemp = Rinhomtemp + ...
                            zonotope([eAt_block * center(project(options.Raux,s_j:f_j)),...
                            eAt_block * generators(project(options.Raux,s_j:f_j))]);
                        Yinhomtemp = Yinhomtemp + ...
                            zonotope([(C{bi} * eAt_block) * center(project(options.Raux,s_j:f_j)),...
                            (C{bi} * eAt_block) * generators(project(options.Raux,s_j:f_j))]);
                    end
                end
            end
            % sum over all bj
            Rinhom{bi} = Rinhomtemp;
            Yinhom{bi} = Yinhomtemp;
        else
            % if no inhomogenuity (no input)
            Rinhom{bi} = 0;
            Yinhom{bi} = 0;
        end
        % -----------------------------------------------------------------
        
    end
end
% -------------------------------------------------------------------------


% save for main loop in reach_decomp --------------------------------------
options.Cblock  = C;
options.Rhom    = Rhom0;
options.Rhom_tp = Rhom_tp0;
options.Cno0    = Cno0;
options.Raux    = Raux;
options.Rtrans  = Rtrans;
options.Yinhom  = Yinhom;
options.Rinhom  = Rinhom;
% -------------------------------------------------------------------------



%------------- END OF CODE --------------