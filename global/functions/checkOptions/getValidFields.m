function validFields = getValidFields(system,type,alg)
% getValidFields - returns options for given system
%
% Syntax:
%    getValidFields(system,type,alg)
%
% Inputs:
%    system      - system type, e.g. linearSys, nonlinParamSys, ...
%    type        - string: 'mand' ... mandatory, 'opt' ... optional,
%                          'def' ... default set, 'all' ... all
%    alg         - algorithm type (only linearSys): see options.alg
%
% Outputs:
%    validFields - cell of strings
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: @contDynamics: checkOptionsReach.m / checkOptionsSimulate.m
%           @hybridAutomaton: checkOptionsReach.m / checkOptionsSimulate.m
%           @parallelHybridAutomaton: checkOptionsReach.m / checkOptionsSimulate.m
%
% References: 
%   -

% Author:       Mark Wetzlinger, Matthias Althoff
% Written:      06-Mar-2019
% Last update:  17-Aug-2019
%               17-Jul-2020 (linProbSys added, MA)
%               24-Jul-2020 (boxInputs removed, MA)
% Last revision:---

%------------- BEGIN CODE --------------

if nargin <= 2
    alg = 'standard';
end
if nargin == 1
    % no 'type'-flag given -> return all options
    type = 'all';
else
    % check if 'type'-flag allowed
    lists = {'mand';'opt';'def';'all'};
    if ~any(strcmp(type,lists))
        error("Incorrent options list definition");
    end
end

% choose correct options lists --------------------------------------------
switch(system)
    case 'contDyn_sim'
        mandList = {'tFinal';'R0';'U';'paramInt';'points'};
        if strcmp(alg,'simRand')
            mandList = [mandList;'fracVert';'fracInpVert';'inpChanges'];
        elseif strcmp(alg,'RRT')
            mandList = [mandList;'vertSamp';'strechFac'];
        end
        optList = {'u';'timeStep'};
        defList = {'tStart'};
        
    case 'linearSys'
        mandList = {'tFinal';'R0';'U'};
        optList = {'u';'saveOrder'};
        defList = {'tStart';'reductionTechnique';'linAlg';'verbose'};
        if ~strcmp(alg,'adap')
            mandList = [mandList; 'timeStep';'taylorTerms';'zonotopeOrder'];
            if strcmp(alg,'decomp')
                mandList = [mandList; 'partition'];
            elseif strcmp(alg,'krylov')
                mandList = [mandList; 'krylovError'; 'krylovStep'];
            end
        else
            defList = [defList; 'error'];
        end
        
    case 'linParamSys'
        mandList = {'tFinal';'R0';'timeStep';'taylorTerms';...
                    'zonotopeOrder';'U';'intermediateOrder'};
        optList = {'u';'saveOrder'};
        defList = {'tStart';'reductionTechnique';'verbose';'compTimePoint'};
        
    case 'linearSysDT'
        mandList = {'tFinal';'R0';'zonotopeOrder';'U'};
        optList = {'u';'saveOrder'};
        defList = {'tStart';'reductionTechnique';'verbose'};
        
    case 'linProbSys'
        mandList = {'tFinal';'R0';'timeStep';'taylorTerms';...
                    'zonotopeOrder';'U';'gamma'};
        optList = {'u';'saveOrder'};
        defList = {'tStart';'reductionTechnique';'linAlg';'verbose';'compTimePoint'};
        
    case 'nonlinDASys'
        mandList = {'tFinal';'R0';'timeStep';'taylorTerms';...
                    'zonotopeOrder';'U';'tensorOrder';...
                    'maxError_x';'maxError_y';'y0guess';'x0'};
        optList = {'saveOrder';'u'
            'intermediateOrder';'errorOrder';'lagrangeRem'};
        defList = {'tStart';'reductionTechnique';'verbose';'linAlg'; ...
                   'maxError';'reductionInterval'};

    case 'nonlinearSysDT'
        mandList = {'tFinal';'R0';'zonotopeOrder'; ...
                    'U';'tensorOrder';'errorOrder'};
        optList = {'lagrangeRem'};
        defList = {'tStart';'reductionTechnique';'verbose'};

    case 'nonlinearSys'
        mandList = {'tFinal';'R0';'timeStep';'taylorTerms';...
                    'zonotopeOrder';'U';'tensorOrder';'alg'};
        optList = {'saveOrder'; 'u'; 'intermediateOrder';'errorOrder'; ...
                   'lagrangeRem'; 'polyZono'};
        defList = {'tStart';'reductionTechnique';'verbose';'linAlg'; ...
                   'maxError';'reductionInterval'};
        
    case 'nonlinParamSys'
        mandList = {'tFinal';'R0';'timeStep';'taylorTerms';...
            'zonotopeOrder';'U';'tensorOrder';'alg';'paramInt'};
        optList = {'saveOrder'; 'u';'intermediateOrder';'errorOrder'; ...
                  'lagrangeRem'; 'polyZono'};
        defList = {'tStart';'reductionTechnique';'verbose';'linAlg'; ...
                   'maxError';'reductionInterval'};

    case 'hybridAutomaton_reach'
        mandList = {'timeStepLoc';'Uloc';'guardIntersect';'startLoc'; ...
                    'finalLoc';'tStart';'tFinal';'enclose'};
        optList = {'intersectInvariant'};
        defList = {};
    
    case 'hybridAutomaton_sim'
        mandList = {'x0';'uLoc';'tFinal';'tStart';'finalLoc';'startLoc'};
        optList = {};
        defList = {};
        
    case 'hybridAutomaton_simRand'
        mandList = {'R0';'Uloc';'tFinal';'tStart';'finalLoc';'startLoc';...
                    'points';'fracVert';'fracInpVert';'inpChanges'};
        optList = {};
        defList = {};
        
    case 'parallelHybridAutomaton_reach'
        mandList = {'timeStep';'tStart';'tFinal';'startLoc';'finalLoc'; ...
                    'Uloc';'inputCompMap';'guardIntersect';'enclose';'R0'};
        optList = {'intersectInvariant'};
        defList = {'linAlg'};
        
    case 'parallelHybridAutomaton_simRand'
        mandList = {'Uloc';'inputCompMap';'startLoc';'finalLoc';'R0'; ...
                    'tStart';'tFinal';...
                    'points';'fracVert';'fracInpVert';'inpChanges'};
        optList = {};
        defList = {};
        
    case 'parallelHybridAutomaton_sim'
        mandList = {'uloc';'inputCompMap';'startLoc';'finalLoc';'x0'; ...
                    'tStart';'tFinal'};
        optList = {};
        defList = {};
        
    otherwise
        error('failed call: incorrect system name');
end
% -------------------------------------------------------------------------


% concatenate lists depending on 'type'-flag ------------------------------
if strcmp(type, 'all')
    validFields = [mandList; optList; defList];
elseif strcmp(type,'mand')
    validFields = mandList;
elseif strcmp(type, 'opt')
    validFields = optList;
elseif strcmp(type, 'def')
    validFields = defList;
end
% -------------------------------------------------------------------------


end


%------------- END OF CODE --------------
