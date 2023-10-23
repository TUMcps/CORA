function defValue = getDefaultValueOptions(field,sys,params,options)
% getDefaultValueOptions - contains list of default values for options
%
% Syntax:
%    defValue = getDefaultValueOptions(field,sys,params,options,listname)
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

% search for default value in options.<field>
switch field
    case 'reductionTechnique'
        defValue = 'girard';
    case 'reductionTechniqueUnderApprox'
        defValue = 'sum';
    case 'linAlg'
        defValue = 'standard';
    case 'zonotopeOrder'
        defValue = inf;
    case 'reachAlg'
        defValue = 'standard';
    case 'verbose'
        defValue = false;
    case 'reductionInterval'
        defValue = Inf;
    case 'maxError'
        defValue = aux_def_maxError(sys,params,options);
    case 'maxError_x'
        defValue = aux_def_maxError_x(sys,params,options);
    case 'maxError_y'
        defValue = aux_def_maxError_y(sys,params,options);
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
        defValue = aux_def_nrConstInp(sys,params,options);
    case 'p_conf'
        defValue = 0.8;
    case 'inpChanges'
        defValue = 1;
    case 'alg'
        defValue = aux_def_alg(sys,params,options);
    case 'tensorOrderOutput'
        defValue = 2;
    case 'compOutputSet'
        defValue = true;
    case 'intersectInvariant'
        defValue = false;
    case 'timeStepDivider'
        defValue = 1;
    case 'postProcessingOrder'
        defValue = inf;
    case 'armaxAlg'
        defValue = 'tvpGeneral';
    case 'norm'
        defValue = 'interval';
    otherwise
        throw(CORAerror('CORA:specialError',...
            "There is no default value for options." + field + "."))
end

end


% Auxiliary functions -----------------------------------------------------

function val = aux_def_maxError(sys,params,options)

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

function val = aux_def_nrConstInp(sys,params,options)

val = [];
if isa(sys,'contDynamics')
    if isa(sys,'linearSysDT') ...
        || isa(sys,'nonlinearSysDT') ...
        || isa(sys,'neurNetContrSys')

        steps = round((params.tFinal - params.tStart) / sys.dt);
        
        if isinf(steps)
            val = 1;
        else
            % start at 10, go down to 1
            for i=10:-1:1
                if mod(steps,i) == 0
                    val = i; break
                end
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

function val = aux_def_maxError_x(sys,params,options)
% only DA systems

val = [];
if isa(sys,'nonlinDASys')
    val = Inf(sys.dim,1);
end

end

function val = aux_def_maxError_y(sys,params,options)
% only DA systems

val = [];
if isa(sys,'nonlinDASys')
    val = Inf(sys.nrOfConstraints,1);
end

end

function val = aux_def_alg(sys,params,options)

% default
val = 'lin';

% explicitly stated for the below classes
if isa(sys,'nonlinearSysDT')
    val = 'lin';
elseif isa(sys,'nonlinDASys')
    val = 'lin';
end

end

% ------------------------------ END OF CODE ------------------------------
