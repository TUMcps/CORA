function defValue = getDefaultValue(field,sys,listname)
% getDefaultValue - contains list of default values for params / options
%
% Syntax:
%    defValue = getDefaultValue(field,sys,listname)
%
% Inputs:
%    field - struct field in params / options
%    sys - object of system class
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

% syntax:
% {'name of param/option #1', default value #1;
%  'name of param/option #2', default value #2;
%    ..., ...};

% define default value for params
params_defFields = {'tStart', 0;
             'finalLoc', def_finalLoc(sys);
             'U', def_U(sys);
             'u', def_u(sys);
             'y', def_y(sys);
             'inputCompMap', def_inputCompMap(sys)};

% define default value for options
options_defFields = {'reductionTechnique', 'girard';
             'reductionTechniqueUnderApprox', 'sum';
             'linAlg', 'standard';
             'verbose', false;
             'reductionInterval', Inf;
             'maxError', def_maxError(sys);
             'compTimePoint', true;
             'polyZono.maxDepGenOrder', 20;
             'polyZono.maxPolyZonoRatio', Inf;
             'polyZono.restructureTechnique', 'reduceGirard';
             'lagrangeRem.simplify', 'none';
             'lagrangeRem.method', 'interval';
             'lagrangeRem.tensorParallel', false;
             'lagrangeRem.optMethod', 'int';
             'linAlg', 'standard';
             'contractor', 'linearize';
             'iter', 2;
             'splits', 8;
             'orderInner', 5;
             'scaleFac', 'auto';
             'inpChanges', 0};
         

% search for default value in params
if strcmp(listname,'params')

    [row,~] = find(strcmp(params_defFields,field));
    if ~isempty(row)
        defValue = params_defFields{row,2}; return
    else
        error("There is no default value for params." + field + "!");
    end
    
% search for default value in params
elseif strcmp(listname,'options')
    
    [row,~] = find(strcmp(options_defFields,field));
    if ~isempty(row)
        defValue = options_defFields{row,2}; return
    else
        error("There is no default value for options." + field + "!");
    end
    
else
    
    error("Not specified which list the field is part of.");
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

end

function val = def_U(sys)

if isa(sys,'contDynamics')
    if isa(sys,'linearSys') || isa(sys,'linearSysDT') || isa(sys,'linParamSys')
        val = zonotope(zeros(size(sys.B,2),1));
    else 
        val = zonotope(zeros(sys.nrOfInputs,1));
    end
elseif isa(sys,'hybridAutomaton')
    locations = sys.location;
    numLoc = length(locations);
    val = cell(numLoc,1);
    for i = 1:numLoc
        subsys = locations{i}.contDynamics;
        val{i} = zonotope(zeros(max(1,subsys.nrOfInputs),1));
    end
elseif isa(sys,'parallelHybridAutomaton')
    numComps = length(sys.components);
    val = cell(numComps,1);
    for i=1:numComps
        numLoc = length(sys.components{i}.location);
        val{i} = cell(numLoc,1);
        for j = 1:numLoc
            val{i}{j} = zonotope(zeros(sys.numInputs,1));
        end    
    end
end

end

function val = def_u(sys)

val = [];
if isa(sys,'contDynamics')
    if isa(sys,'linearSys') || isa(sys,'linearSysDT') || isa(sys,'linParamSys') || isa(sys,'linProbSys')
        val = zeros(size(sys.B,2),1);
    else 
        val = zeros(sys.nrOfInputs,1);
    end
elseif isa(sys,'hybridAutomaton')
    locations = sys.location;
    numLoc = length(locations);
    val = cell(numLoc,1);
    for i = 1:numLoc
        subsys = locations{i}.contDynamics;
        val{i} = zeros(max(1,subsys.nrOfInputs),1);
    end
elseif isa(sys,'parallelHybridAutomaton')
    numComps = length(sys.components);
    val = cell(numComps,1);
    for i=1:numComps
        numLoc = length(sys.components{i}.location);
        val{i} = cell(numLoc,1);
        for j = 1:numLoc
            val{i}{j} = zeros(sys.numInputs,1);
        end    
    end
end

end

function val = def_y(sys)

val = [];
if isa(sys,'contDynamics')
    val = zeros(sys.nrOfOutputs,1);
end
% no assignment for hybridAutomaton / parallelHybridAutomaton

end

function val = def_inputCompMap(sys)

val = [];
if isa(sys,'parallelHybridAutomaton')
    val = ones(sys.numInputs,1);
end

end

function val = def_maxError(sys)

val = [];
if isa(sys,'contDynamics')
    val = Inf(sys.dim,1);
end

end


%------------- END OF CODE --------------