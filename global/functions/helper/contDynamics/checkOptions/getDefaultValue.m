function defValue = getDefaultValue(field,sys,params,options,listname)
% getDefaultValue - contains list of default values for params / options
%
% Syntax:
%    defValue = getDefaultValue(field,sys,params,options,listname)
%
% Inputs:
%    field - struct field in params / options
%    sys - object of system class
%    params - struct containing model parameters
%    options - struct containing algorithm parameters
%    listname - 'params' or 'options'
%
% Outputs:
%    defValue - default value for given field
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      26-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% split list for params / options for more transparency / readability
switch listname
    
    % search for default value in params
    case 'params'
    
    switch field
        case 'tStart'
            defValue = 0;
        case 'finalLoc'
            defValue = def_finalLoc(sys);
        case 'U'
            defValue = def_U(sys);
        case 'u'
            defValue = def_u(sys);
        case 'tu'
            defValue = def_tu(sys,params);
        case 'y'
            defValue = def_y(sys);
        case 'W'
            defValue = def_W(sys);
        case 'V'
            defValue = def_V(sys);
        case 'inputCompMap'
            defValue = def_inputCompMap(sys);
        otherwise
            throw(CORAerror('CORA:specialError',...
                "There is no default value for params." + field + "."))
    end
    
    % search for default value in options
    case 'options'
        
        switch field
            case 'reductionTechnique'
                defValue = 'girard';
            case 'reductionTechniqueUnderApprox'
                defValue = 'sum';
            case 'linAlg'
                defValue = 'standard';
            case 'verbose'
                defValue = false;
            case 'reductionInterval'
                defValue = Inf;
            case 'maxError'
                defValue = def_maxError(sys);
            case 'maxError_x'
                defValue = def_maxError_x(sys);
            case 'maxError_y'
                defValue = def_maxError_y(sys);
            case 'compTimePoint'
                defValue = true;
            case 'polyZono.maxDepGenOrder'
                defValue = 20;
            case 'polyZono.maxPolyZonoRatio'
                defValue = Inf;
            case 'polyZono.restructureTechnique'
                defValue = 'reduceGirard';
            case 'lagrangeRem.simplify'
                defValue = 'none';
            case 'lagrangeRem.method'
                defValue = 'interval';
            case 'lagrangeRem.tensorParallel'
                defValue = false;
            case 'lagrangeRem.optMethod'
                defValue = 'int';
            case 'contractor'
                defValue = 'linearize';
            case 'iter'
                defValue = 2;
            case 'splits'
                defValue = 8;
            case 'orderInner'
                defValue = 5;
            case 'scaleFac'
                defValue = 'auto';
            case 'type'
                defValue = 'standard';
            case 'points'
                defValue = 10;
            case 'fracVert'
                defValue = 0.5;
            case 'fracInpVert'
                defValue = 0.5;
            case 'nrConstInp'
                defValue = def_nrConstInp(sys,params);
            case 'p_conf'
                defValue = 0.8;
            case 'inpChanges'
                defValue = 1;
            case 'alg'
                defValue = def_alg(sys);
            case 'tensorOrderOutput'
                defValue = 2;
            case 'compOutputSet'
                defValue = true;
            case 'intersectInvariant'
                defValue = false;
            otherwise
                throw(CORAerror('CORA:specialError',...
                    "There is no default value for options." + field + "."))
        end
end

end


% Auxiliary Functions for split between hybrid and contDynamics classes
function val = def_finalLoc(sys)

val = [];
if isa(sys,'hybridAutomaton')
    val = 0;
elseif isa(sys,'parallelHybridAutomaton')
    val = zeros(length(sys.components),1);
end
% no assignment for contDynamics

end

function val = def_U(sys)

if isa(sys,'contDynamics')
    val = zonotope(zeros(sys.nrOfInputs,1));
elseif isa(sys,'hybridAutomaton')
    locations = sys.location;
    numLoc = length(locations);
    val = cell(numLoc,1);
    for i = 1:numLoc
        nrInputs = locations(i).contDynamics.nrOfInputs;
        val{i} = zonotope(zeros(max(1,nrInputs),1));
    end
elseif isa(sys,'parallelHybridAutomaton')
    numComps = length(sys.components);
    val = cell(numComps,1);
    for i=1:numComps
        numLoc = length(sys.components(i).location);
        val{i} = cell(numLoc,1);
        for j = 1:numLoc
            nrInputs = sys.components(i).location(1).contDynamics.nrOfInputs;
            val{i}{j} = zonotope(zeros(nrInputs,1));
        end    
    end
end

end

function val = def_u(sys)

val = [];
if isa(sys,'contDynamics')
    val = zeros(sys.nrOfInputs,1);
elseif isa(sys,'hybridAutomaton')
    locations = sys.location;
    numLoc = length(locations);
    val = cell(numLoc,1);
    for i = 1:numLoc
        subsys = locations(i).contDynamics;
        val{i} = zeros(max(1,subsys.nrOfInputs),1);
    end
elseif isa(sys,'parallelHybridAutomaton')
    numComps = length(sys.components);
    val = cell(numComps,1);
    for i=1:numComps
        numLoc = length(sys.components(i).location);
        val{i} = cell(numLoc,1);
        for j = 1:numLoc
            val{i}{j} = zeros(sys.nrOfInputs,1);
        end    
    end
end

end

function val = def_tu(sys,params)

val = [];
if isa(sys,'contDynamics')
    if isa(sys,'linearSysDT') || isa(sys,'nonlinearSysDT')
        if size(params.u,2) == 1
            val = params.tStart;
        else
            val = (params.tStart:sys.dt:params.tFinal-sys.dt)';
            if isa(sys,'linearSysDT') && any(any(sys.D))
                val = (params.tStart:sys.dt:params.tFinal)';
            end
        end
    elseif isa(sys,'linearSys')
        steps = size(params.u,2);
        if any(any(sys.D)) && steps > 1
            steps = steps - 1;
        end
        stepsize = (params.tFinal-params.tStart) / steps;
        val = (params.tStart:stepsize:params.tFinal-stepsize)';
        if steps > 1 && any(any(sys.D))
            val = (params.tStart:stepsize:params.tFinal)';
        end
    else % isa(sys,'linParamSys') || isa(sys,'linProbSys')
        if size(params.u,2) == 1
            val = params.tStart;
        else
            steps = size(params.u,2);
            stepsize = (params.tFinal-params.tStart) / steps;
%             val = linspace(params.tStart,params.tFinal-stepsize,steps);
            val = (params.tStart:stepsize:params.tFinal-stepsize)';
        end
    end
        
end
% no assignment for hybridAutomaton / parallelHybridAutomaton

end

function val = def_y(sys)

val = [];
if isa(sys,'contDynamics')
    val = zeros(sys.nrOfOutputs,1);
end
% no assignment for hybridAutomaton / parallelHybridAutomaton

end

function val = def_W(sys)

if isa(sys,'contDynamics') || isa(sys,'parallelHybridAutomaton')
    % note: dimension of composed automaton is fixed
    val = interval(zeros(sys.dim,1));
elseif isa(sys,'hybridAutomaton')
    n = sys.dim;
    if all(n(1) == n)
        val = interval(zeros(n(1),1),zeros(n(1),1));
    else
        throw(CORAerror('CORA:notSupported',...
            'Default value for W not supported for hybrid automata with varying number of states per location.'));
    end
end

end

function val = def_V(sys)

if isa(sys,'contDynamics') || isa(sys,'parallelHybridAutomaton')
    % note: dimension of composed automaton is fixed
% CORA does not support output equations for composed parallel hybrid
% automata -> use state dimension instead
    val = interval(zeros(sys.nrOfOutputs,1),zeros(sys.nrOfOutputs,1));
elseif isa(sys,'hybridAutomaton')
    r = sys.nrOfOutputs;
    if all(r(1) == r)
        val = interval(zeros(r(1),1));
    else
        throw(CORAerror('CORA:notSupported',...
            'Default value for W not supported for hybrid automata with varying number of states per location.'));
    end
end

end

function val = def_inputCompMap(sys)

val = [];
if isa(sys,'parallelHybridAutomaton')
    val = ones(sys.nrOfInputs,1);
end
% no assignment for contDynamics

end

function val = def_maxError(sys)

val = [];
if isa(sys,'contDynamics') || isa(sys,'parallelHybridAutomaton')
    val = Inf(sys.dim,1);
elseif isa(sys,'hybridAutomaton')
    n = sys.dim;
    if all(n(1) == n)
        val = Inf(n(1),1);
    else
        throw(CORAerror('CORA:notSupported',...
            'Default value for maxError not supported for hybrid automata with varying number of states per location.'));
    end
end
% no assignment for hybridAutomaton / parallelHybridAutomaton

end

function val = def_nrConstInp(sys,params)

val = [];
if isa(sys,'contDynamics')
    if isa(sys,'linearSysDT') || isa(sys,'nonlinearSysDT') || ...
                                                isa(sys,'neurNetContrSys')
        steps = round((params.tFinal - params.tStart) / sys.dt);
        % start at 10, go down to 1
        for i=10:-1:1
            if mod(steps,i) == 0
                val = i; break
            end
        end
    else
        if size(params.u,2) > 1
            if isa(sys,'linearSys') && any(any(sys.D))
                val = size(params.u,2) - 1;
            else
                val = size(params.u,2);
            end
        else % no input trajectory
            val = 10;
        end
    end
elseif isa(sys,'hybridAutomaton') || isa(sys,'parallelHybridAutomaton')
    % this will most likely be changed in the future (e.g., different
    % values for each location)
    val = 10;
end

end

function val = def_maxError_x(sys)
% only DA systems

val = [];
if isa(sys,'nonlinDASys')
    val = Inf(sys.dim,1);
end

end

function val = def_maxError_y(sys)
% only DA systems

val = [];
if isa(sys,'nonlinDASys')
    val = Inf(sys.nrOfConstraints,1);
end

end

function val = def_alg(sys)

% default
val = 'lin';

% explicitly stated for the below classes
if isa(sys,'nonlinearSysDT')
    val = 'lin';
elseif isa(sys,'nonlinDASys')
    val = 'lin';
end

end

%------------- END OF CODE --------------