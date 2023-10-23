function defValue = getDefaultValueParams(field,sys,params,options)
% getDefaultValueParams - contains list of default values for params
%
% Syntax:
%    defValue = getDefaultValueParams(field,sys,params,options,listname)
%
% Inputs:
%    field - struct field in params / options
%    sys - object of system class
%    params - struct containing model parameters
%    options - struct containing algorithm parameters
%
% Outputs:
%    defValue - default value for given field
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getDefaultValues

% Authors:       Mark Wetzlinger
% Written:       26-January-2021
% Last update:   09-October-2023 (TL, split options/params)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% search for default value in params.<field>
switch field
    case 'tStart'
        defValue = 0;
    case 'finalLoc'
        defValue = aux_def_finalLoc(sys,params,options);
    case 'U'
        defValue = aux_def_U(sys,params,options);
    case 'u'
        defValue = aux_def_u(sys,params,options);
    case 'tu'
        defValue = aux_def_tu(sys,params,options);
    case 'y'
        defValue = aux_def_y(sys,params,options);
    case 'W'
        defValue = aux_def_W(sys,params,options);
    case 'V'
        defValue = aux_def_V(sys,params,options);
    case 'inputCompMap'
        defValue = aux_def_inputCompMap(sys,params,options);
    % conformance
    case 'w'
        defValue = aux_def_conformWeight(sys,params);
    otherwise
        throw(CORAerror('CORA:specialError',...
            "There is no default value for params." + field + "."))
end 

end


% Auxiliary functions -----------------------------------------------------

function val = aux_def_finalLoc(sys,params,options)

val = [];
if isa(sys,'hybridAutomaton')
    val = 0;
elseif isa(sys,'parallelHybridAutomaton')
    val = zeros(length(sys.components),1);
end
% no assignment for contDynamics

end

function val = aux_def_U(sys,params,options)

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

function val = aux_def_u(sys,params,options)

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

function val = aux_def_tu(sys,params,options)

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

function val = aux_def_y(sys,params,options)

val = [];
if isa(sys,'contDynamics')
    val = zeros(sys.nrOfOutputs,1);
end
% no assignment for hybridAutomaton / parallelHybridAutomaton

end

function val = aux_def_W(sys,params,options)

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

function val = aux_def_V(sys,params,options)

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

function val = aux_def_inputCompMap(sys,params,options)

val = [];
if isa(sys,'parallelHybridAutomaton')
    val = ones(sys.nrOfInputs,1);
end
% no assignment for contDynamics

end

function val = aux_def_conformWeight(sys,params,options)
    maxNrOfTimeSteps = ceil(params.tFinal/sys.dt); % maximum number of timeSteps
    val = ones(maxNrOfTimeSteps+1,1);
end

% ------------------------------ END OF CODE ------------------------------
