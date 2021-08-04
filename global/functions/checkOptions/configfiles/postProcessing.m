function [params,options] = postProcessing(sys,func,params,options)
% postProcessing - perform post-processing operations for
%    model parameters and algorithm parameters
%
% Syntax:
%    [params,options] = postProcessing(sys,func,params,options)
%
% Inputs:
%    sys - contDynamics or hybridDynamics object
%    func - function name ('reach', 'simulateRandom', etc.)
%    params - user-defined model parameters
%    options - user-defined algorithm parameters
%
% Outputs:
%    params - user-defined model parameters updated for internal usage
%    options - user-defined algorithm parameters updated for internal usage
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger, Matthias Althoff
% Written:      04-February-2021
% Last update:  07-July-2021 (MA, U should not be converted for simulateGaussian and observe)
% Last revision:---

%------------- BEGIN CODE --------------


% only contDynamics -------------------------------------------------------

if isa(sys,'contDynamics')
    
    % convert U to a zonotope if given as an interval
    if ~any(strcmp(func,{'reachInnerProjection','simulateGaussian','observe'}))
        [params,options] = convert_U(sys,params,options);
    end
    
    if ~isa(sys,'linearSys')
        [params,options] = set_linAlg(sys,params,options);
    end
    
    % extend vectors to correct length if default setting
    if strcmp(func,'observe')
        [params,options] = set_u_y(sys,params,options);
    end
    
    % set input set U and input vector uTrans|uTransVec
    if ~(contains(func,'reachInner') || strcmp(func,'observe'))
        [params,options] = set_U_uTrans_uTransVec(sys,params,options);
    
        % determine originContained
        if ~any(strcmp(func,{'simulateNormal','simulateRRT','simulateGaussian'})) && ...
            (isa(sys,'linearSys') || isa(sys,'linParamSys') || isa(sys,'linProbSys'))
            [params,options] = set_originContained(sys,params,options);
        end
        
    end
    
    % set alg and N
    if isa(sys,'nonlinearSys') && contains(func,'reachInner')
        if strcmp(func,'reachInnerScaling')
            % use polynomialization, compute number of steps
            [params,options] = set_R0_alg_tensorOrder(sys,params,options);
            [params,options] = set_N(sys,params,options);
        elseif strcmp(func,'reachInnerProjection')
            % R0 converted to an interval
            [params,options] = set_R0int(sys,params,options);
        end
    end
    
    % set polyZonotope restructuring
    if ~any(strcmp(func,{'simulateNormal','simulateRRT','simulateGaussian'}))
        if isa(sys,'nonlinearSysDT')
            [params,options] = set_volApproxMethod(sys,params,options,1);
        elseif isa(sys,'nonlinearSys') || isa(sys,'nonlinParamSys')
            [params,options] = set_volApproxMethod(sys,params,options,0);
        end
    end
    
    % nonlinear adaptive tuning
    if strcmp(func,'reach') && isa(sys,'nonlinearSys') && contains(options.alg,'adaptive')
        [params,options] = set_nonlinearAdaptive(sys,params,options);
    end
    
end

% only hybridAutomaton ----------------------------------------------------

if isa(sys,'hybridAutomaton')
   
    if strcmp(func,'reach')
        % set timeStep for each location
        [params,options] = set_HA_timeStep(sys,params,options);
    end
    if any(strcmp(func,{'reach','simulateRandom'}))
        % set input for locations
        [params,options] = set_HA_inputs(sys,params,options);
    end
    if strcmp(func,'simulate')
        % set input for locations
        [params,options] = set_HA_sim_inputs(sys,params,options);
    end
    
end

% only parallelHybridAutomaton --------------------------------------------

if isa(sys,'parallelHybridAutomaton')
   
    % set input for locations
    if any(strcmp(func,{'reach','simulateRandom'}))
        [params,options] = set_pHA_inputs(sys,params,options);
    end
    if strcmp(func,'simulate')
        [params,options] = set_pHA_sim_inputs(sys,params,options);
    end
    
end

end

% Auxiliary Functions: contDynamics

function [params,options] = convert_U(sys,params,options)

% convert U to a zonotope if given as an interval
if isa(params.U,'interval')
    params.U = zonotope(params.U);
end

end

function [params,options] = set_linAlg(sys,params,options)

% only 'standard' for all linearized systems
options.linAlg = 'standard';

end

function [params,options] = set_U_uTrans_uTransVec(sys,params,options)

% shift u by center of U
centerU = center(params.U);
if any(centerU)
    params.u = params.u + center(params.U);
    params.U = params.U + (-centerU);
end

% write u to uTrans or uTrans
if size(params.u,2) > 1
    params.uTransVec = params.u;
else
    params.uTrans = params.u;
end
params = rmfield(params,'u');

end

function [params,options] = set_originContained(sys,params,options)

U = params.U;
% define actual input set
if isfield(params,'uTrans')
    U = U + params.uTrans;
elseif isfield(params,'uTransVec')
    % input trajectory
    options.originContained = false; return;
end

% options.originContained = false;
if isa(sys,'linearSys') && ~isempty(sys.c)
    % check if any of equality constraints B(i)*U = -c(i) = 0
    % is unsatisfiable
    found = false;
    for i = 1:size(sys.B,1)
       hp = conHyperplane(sys.B(i,:),-sys.c(i));
       if ~isIntersecting(hp,U)
          options.originContained = false;
          found = true;
       end
    end
    if ~found
        % check if vTrans = B*U + c contains the origin
        vTrans = sys.B*U + sys.c;
        % faster computation if vTrans is an interval
        if isInterval(vTrans)
            vTransInt = interval(vTrans);
            options.originContained = in(vTransInt,zeros(dim(vTrans),1));
        else
            options.originContained = in(vTrans,zeros(dim(vTrans),1));
        end
    end
elseif isInterval(U)
    % no constant input, faster computation if U is an interval
    int = interval(U);
    options.originContained = in(int,zeros(dim(U),1));
else
    % no constant input c, U not an interval
    options.originContained = in(U,zeros(dim(U),1));
end

end

function [params,options] = set_u_y(sys,params,options)

if size(params.u,2) == 1 && ~any(params.u)
    % extend vector to correct number of columns
    params.u = repmat(params.u,1,(params.tFinal-params.tStart)/sys.dt);   
end
% rewrite to uTransVec
params.uTransVec = params.u;
params = rmfield(params,'u');

if size(params.y,2) == 1 && ~any(params.y)
    % extend vector to correct number of columns
    params.y = repmat(params.y,1,(params.tFinal-params.tStart)/sys.dt);   
end

end

function [params,options] = set_R0_alg_tensorOrder(sys,params,options)

params.R0 = polyZonotope(params.R0);
options.alg = 'poly';
options.tensorOrder = 3;

end

function [params,options] = set_R0int(sys,params,options)

params.R0 = interval(params.R0);

end

function [params,options] = set_volApproxMethod(sys,params,options,isDT)

if ~isDT
    % isfield because of nonlinearSys/reachInner
    if isfield(options,'alg') && strcmp(options.alg,'poly')
        if isFullDim(params.R0)
            options.polyZono.volApproxMethod = 'interval';
        else
            options.polyZono.volApproxMethod = 'pca';
        end
    end
else
    if isFullDim(params.R0)
        options.polyZono.volApproxMethod = 'interval';
    else
        options.polyZono.volApproxMethod = 'pca';
    end
end

end

function [params,options] = set_N(sys,params,options)

if isfield(options,'timeStepInner')
    options.N = round(options.timeStepInner/options.timeStep,0);
    options = rmfield(options,'timeStepInner');
else
    options.N = length(params.tStart:options.timeStep:params.tFinal)-1;
end

end

function [params,options] = set_nonlinearAdaptive(sys,params,options)

options.decrFactor = 0.90;              % adaptation of delta t
options.zetaTlin = 0.0005;              % zeta_T,lin (taylorTerms)
options.zetaTabs = 0.005;               % zeta_T,abs (taylorTerms)
if contains(options.alg,'lin')
    params.R0 = zonotope(params.R0);    % convert to zono
    params.R = params.R0;               % for initial set
    options.redFactor = 0.0005;         % zeta_Z (zonotope order)
    options.zetaK = 0.90;               % zeta_K (tensorOrder)
    options.zetaphi = 0.85;             % zeta_Delta (timeStep)
elseif contains(options.alg,'poly')
    % options for poly (fixed in adaptive)
    options.volApproxMethod = 'interval';           % polyZono
    options.maxDepGenOrder = 50;                    % polyZono
    options.maxPolyZonoRatio = 0.05;                % polyZono
    options.restructureTechnique = 'reducePca';     % polyZono
    params.R0 = polyZonotope(params.R0);            % convert to polyZono
    params.R = params.R0;                           % for initial set
    options.redFactor = 0.0002;                     % zeta_Z (pZ order)
    options.zetaphi = 0.80;                         % zeta_Delta (timeStep)
    options.tensorOrder = 3;                        % fixed
end
options.R.error = zeros(sys.dim,1);                 % for consistency

% options to speed up tensor computation
options.thirdOrderTensorempty = false;
options.isHessianConst = false;
options.hessianCheck = false;

end

% Auxiliary Functions: hybridAutomaton and parallelHybridAutomaton

function [params,options] = set_HA_timeStep(sys,params,options)
% timeStep becomes timeStepLoc for each location ... only if non-adaptive

if ~strcmp(options.linAlg,'adaptive')
    
    if ~iscell(options.timeStep)
        locations = sys.location;
        numLoc = length(locations);
        options.timeStepLoc = cell(numLoc,1);
        for i=1:numLoc
            options.timeStepLoc{i} = options.timeStep;
        end
    else
        options.timeStepLoc = options.timeStep;
    end
    
    options = rmfield(options,'timeStep');
    
end

end

function [params,options] = set_HA_inputs(sys,params,options)
% internal value options.Uloc stores the input set for each location

locations = sys.location;
numLoc = length(locations);

if ~iscell(params.U)
    isaInt = false;
    if isa(params.U,'interval')
        isaInt = true; UZon = zonotope(params.U);
    end
    % same input set for each location
    Uloc = cell(numLoc,1);
    for i = 1:numLoc
        if isaInt
            Uloc{i} = UZon;
        else
            Uloc{i} = params.U;
        end
    end
    
else
    % copy input set, convert to zonotope if necessary
    Uloc = cell(numLoc,1);
    for i = 1:numLoc
        if isa(params.U{i},'interval')
            Uloc{i} = zonotope(params.U{i});
        else
            Uloc{i} = params.U{i};
        end
    end
end

params = rmfield(params,'U');
params.Uloc = Uloc;

end

function [params,options] = set_HA_sim_inputs(sys,params,options)
% internal value options.uLoc which the input set for each location

locations = sys.location;
numLoc = length(locations);

if ~iscell(params.u) 
    % same input set for each location    
    uloc = cell(numLoc,1);
    for i = 1:numLoc
        uloc{i} = params.u;
    end
else
    % copy input set for each location
    uloc = params.u;
end

params = rmfield(params,'u');
params.uLoc = uloc;

end

function [params,options] = set_pHA_inputs(sys,params,options)
% internal value options.Uloc stores the input set for each location

numComps = length(sys.components);

if ~iscell(params.U)
    isaInt = false;
    if isa(params.U,'interval')
        isaInt = true; UZon = zonotope(params.U);
    end
    % same input set for each location
    uloc = cell(numComps,1);
    numLoc = length(sys.components{1}.location);
    uloc{1} = cell(numLoc,1);    
    for i = 1:numLoc
        if isaInt
            uloc{1}{i} = UZon;
        else
            uloc{1}{i} = params.U;
        end
    end    
    params.Uloc = uloc;
else
    % convert to zonotope 
    for i = 1:numComps
        numLoc = length(sys.components{1}.location);
        for j = 1:numLoc
            if isa(params.U{i}{j},'interval')
                uloc{i}{j} = zonotope(params.U{i}{j});
            else
                uloc{i}{j} = params.U{i}{j};
            end
        end
    end    
    params.Uloc = uloc;
end

params = rmfield(params,'U');

end

function [params,options] = set_pHA_sim_inputs(sys,params,options)
% internal value options.uloc stores the input set for each location

numComps = length(sys.components);

if ~iscell(params.u) 
    % same input set for each location
    uloc = cell(numComps,1);
    numLoc = length(sys.components{1}.location);
    uloc{1} = cell(numLoc,1);
    for i = 1:numLoc
        uloc{1}{i} = params.u;
    end   
    params.uLoc = uloc;
else
    params.uLoc = params.u;
end

params = rmfield(params,'u');

end



%------------- END OF CODE --------------

