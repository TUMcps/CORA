function checks = checkDynParameterParams(field,sys,func,params,options,checks)
% checkDynParameterParams - checks dynamic parameter values
%
% Syntax:
%    checkDynParameterParams(field,sys,params,options,checks)
%
% Inputs:
%    field - struct field in params / options
%    sys - object of system class
%    func - function
%    params - struct containing model parameters
%    options - struct containing algorithm parameters
%    checks - struct
%
% Outputs:
%    checks - struct
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: checkDynParameter

% Authors:       Tobias Ladner
% Written:       05-October-2023
% Last update:   09-October-2023 (TL, split options/params)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% search for checks in params
switch field
    case 'tStart'
        checks = aux_getChecksParams_tStart(checks,sys,func,params,options);
    case 'tFinal'
        checks = aux_getChecksParams_tFinal(checks,sys,func,params,options);
    case 'R0'
        checks = aux_getChecksParams_R0(checks,sys,func,params,options);
    case 'U'
        checks = aux_getChecksParams_U(checks,sys,func,params,options);
    case 'u'
        checks = aux_getChecksParams_u(checks,sys,func,params,options);
    case 'tu'
        checks = aux_getChecksParams_tu(checks,sys,func,params,options);
    case 'W'
        checks = aux_getChecksParams_W(checks,sys,func,params,options);
    case 'V'
        checks = aux_getChecksParams_V(checks,sys,func,params,options);
    case 'y'
        checks = aux_getChecksParams_y(checks,sys,func,params,options);
    case 'safeSet'
        checks = aux_getChecksParams_safeSet(checks,sys,func,params,options);
    case 'unsafeSet'
        checks = aux_getChecksParams_unsafeSet(checks,sys,func,params,options);
    case 'paramInt'
        checks = aux_getChecksParams_paramInt(checks,sys,func,params,options);
    case 'y0guess'
        checks = aux_getChecksParams_y0guess(checks,sys,func,params,options);
    case 'startLoc'
        checks = aux_getChecksParams_startLoc(checks,sys,func,params,options);
    case 'finalLoc'
        checks = aux_getChecksParams_finalLoc(checks,sys,func,params,options);
    case 'x0'
        checks = aux_getChecksParams_x0(checks,sys,func,params,options);
    case 'refPoints'
        checks = aux_getChecksParams_refPoints(checks,sys,func,params,options);
    case 'paramInts'
        checks = aux_getChecksParams_paramInts(checks,sys,func,params,options);
    case 'inputCompMap'
        checks = aux_getChecksParams_inputCompMap(checks,sys,func,params,options);
    case 'R0conf'
        checks = aux_getChecksParams_R0conf(checks,sys,func,params,options);
    case 'testSuite'
        checks = aux_getChecksParams_testSuite(checks,sys,func,params,options);
    case 'w'
        checks = aux_getChecksParams_w(checks,sys,func,params,options);
    case 'y0'
        checks = aux_getChecksParams_y0(checks,sys,func,params,options);

    otherwise
        warning('CORA: Unknown params.%s', field); return;
end

end


% Auxiliary functions -----------------------------------------------------

% params.<field> ----------------------------------------------------------

% tStart
function checks = aux_getChecksParams_tStart(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@(val)ge(val,0), 'gezero');
end

% tFinal
function checks = aux_getChecksParams_tFinal(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isscalar, 'isscalar');
    checks(end+1) = add2checks(@(val)ge(val,params.tStart), 'getStart');
end

% R0
function checks = aux_getChecksParams_R0(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('R0'),class(val))), 'memberR0');
    if isa(sys,'hybridAutomaton')
        checks(end+1) = add2checks(@(val)dim(val)==sys.location(params.startLoc).contDynamics.dim, 'eqsysdim');
    else
        checks(end+1) = add2checks(@(val)eq(dim(val),sys.dim), 'eqsysdim');
    end
end

% U
function checks = aux_getChecksParams_U(checks,sys,func,params,options)
    if isa(sys,'hybridAutomaton')
        checks(end+1) = add2checks(@(val)(~iscell(val)&&isa(val,'contSet'))||(iscell(val)&&(all(size(val)==[length(sys.location),1])||all(size(val)==[1,length(sys.location)]))), '???');
    elseif isa(sys,'parallelHybridAutomaton')
        checks(end+1) = add2checks(@(val)c_pHA_U(val,sys,params),'');
    else
        checks(end+1) = add2checks(@(val)any(ismember(getMembers('U'),class(val))), 'memberU');
        checks(end+1) = add2checks(@(val)eq(dim(val),sys.nrOfInputs), 'eqinput');
    end
end

% u
function checks = aux_getChecksParams_u(checks,sys,func,params,options)
    if isa(sys,'hybridAutomaton')
        if strcmp(func,'simulate')
            checks(end+1) = add2checks(@(val)c_HA_sim_u(val,sys,params), '');
        else
            checks(end+1) = add2checks(@(val)all(size(val) == [length(sys.location),1]), '???');
        end
    elseif isa(sys,'parallelHybridAutomaton')
        if strcmp(func,'simulate')
            checks(end+1) = add2checks(@(val)c_pHA_sim_u(val,sys,params), '');
        else
            checks(end+1) = add2checks(@iscell, 'iscell');
            checks(end+1) = add2checks(@(val)all(size(val) == [length(sys.components),1]), '???');
        end
    else
        checks(end+1) = add2checks(@isnumeric, 'isnumeric');
        checks(end+1) = add2checks(@(val)eq(size(val,1),sys.nrOfInputs), 'eqinput');
    end
end

% tu
function checks = aux_getChecksParams_tu(checks,sys,func,params,options)
    if isa(sys,'hybridAutomaton')
        % no check
    else
        checks(end+1) = add2checks(@isvector, 'isvector');
        checks(end+1) = add2checks(@isnumeric, 'isnumeric');
        checks(end+1) = add2checks(@(val)all(diff(val)>0), 'vectorgezero');
        checks(end+1) = add2checks(@(val)length(val)==size(params.u,2), 'equ');
        checks(end+1) = add2checks(@(val)c_tu(val,sys,params,options), '');
    end
end

% W
function checks = aux_getChecksParams_W(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('W'),class(val))), 'memberW');
    checks(end+1) = add2checks(@(val)eq(dim(val),sys.dim), 'eqsysdim');
end

% V
function checks = aux_getChecksParams_V(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('V'),class(val))), 'memberW');
    checks(end+1) = add2checks(@(val)eq(dim(val),sys.nrOfOutputs), 'eqoutput');
end

% y
function checks = aux_getChecksParams_y(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)c_measurements(val,sys,params), '');
end

% safeSet
function checks = aux_getChecksParams_safeSet(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)c_safeSet(val,sys,params,options), '');
end

% unsafeSet
function checks = aux_getChecksParams_unsafeSet(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)c_unsafeSet(val,sys,params,options), '');
end

% paramInt
function checks = aux_getChecksParams_paramInt(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)length(val)==sys.nrOfParam, 'eqparam');
    checks(end+1) = add2checks( @(val)isa(val,'interval') || (isvector(val) && isnumeric(val)), 'vectororinterval');
end

% y0guess
function checks = aux_getChecksParams_y0guess(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)length(val)==sys.nrOfConstraints, 'eqconstr');
end

% startLoc
function checks = aux_getChecksParams_startLoc(checks,sys,func,params,options)
    if isa(sys,'hybridAutomaton')
        checks(end+1) = add2checks(@isscalar, 'isscalar');
        checks(end+1) = add2checks(@isnumeric, 'isnumeric');
        checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
        checks(end+1) = add2checks(@(val)ge(val,0), 'geone');
        checks(end+1) = add2checks(@(val)le(val,length(sys.location)), 'leloc');
    else
        checks(end+1) = add2checks(@(val)c_pHA_startLoc(val,sys,params), '');
    end
end

% finalLoc
function checks = aux_getChecksParams_finalLoc(checks,sys,func,params,options)
    if isa(sys,'hybridAutomaton')
        checks(end+1) = add2checks(@isscalar, 'isscalar');
        checks(end+1) = add2checks(@isnumeric, 'isnumeric');
        checks(end+1) = add2checks(@(val)mod(val,1)==0, 'integer');
        checks(end+1) = add2checks(@(val)ge(val,0), 'geone');
        checks(end+1) = add2checks(@(val)le(val,length(sys.location)+1), 'lelocplus1');
    else
        checks(end+1) = add2checks(@(val)c_pHA_finalLoc(val,sys,params), '');
    end
end

% x0
function checks = aux_getChecksParams_x0(checks,sys,func,params,options)
    if isa(sys,'hybridAutomaton')
        checks(end+1) = add2checks(@isvector, 'isvector');
        checks(end+1) = add2checks(@(val)length(val)==sys.location(params.startLoc).contDynamics.dim, 'isnumeric');
    elseif isa(sys,'parallelHybridAutomaton')
        checks(end+1) = add2checks(@(val)all(size(val)==[sys.dim,1]), 'eqsysdim');
    end
end

% refPoints
function checks = aux_getChecksParams_refPoints(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)size(val,1)==sys.dim, 'eqsysdim');
    checks(end+1) = add2checks(@(val)size(val,2)==reachSteps(params,options)+1, 'eqreachSteps');
end

% paramInts
function checks = aux_getChecksParams_paramInts(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)size(val,1)==sys.dim, 'eqsysdim');
    checks(end+1) = add2checks(@(val)size(val,2)==reachSteps(params,options)+1, 'eqreachSteps');
end

% inputCompMap
function checks = aux_getChecksParams_inputCompMap(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)max(val)<=length(sys.components), 'lecomp');
    checks(end+1) = add2checks(@(val)min(val)>=1, 'vectorgeone');
    checks(end+1) = add2checks(@(val) all(size(val)==[sys.nrOfInputs,1]) || all(size(val)==[1,sys.nrOfInputs]), '???');
end

% R0conf
function checks = aux_getChecksParams_R0conf(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)any(ismember(getMembers('R0conf'),class(val))), 'memberR0conf');
    checks(end+1) = add2checks(@(val)eq(dim(val),sys.dim), 'eqsysdim');
end

% testSuite
function checks = aux_getChecksParams_testSuite(checks,sys,func,params,options)
    checks(end+1) = add2checks(@(val)all(cellfun(@(x)isa(x,'testCase'),val)), 'istestCase');
end

% w
function checks = aux_getChecksParams_w(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isvector, 'isvector');
    checks(end+1) = add2checks(@(x) all(x>=0), 'vectorgezero');
end

% y0
function checks = aux_getChecksParams_y0(checks,sys,func,params,options)
    checks(end+1) = add2checks(@isnumeric, 'isnumeric');
    checks(end+1) = add2checks(@(val)eq(size(val,1),sys.nrOfOutputs), 'eqoutput');
end

% ------------------------------ END OF CODE ------------------------------
