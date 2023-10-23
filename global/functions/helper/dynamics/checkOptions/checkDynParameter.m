function failedChecks = checkDynParameter(field,sys,func,params,options,listname)
% checkDynParameter - checks dynamic parameter values
%
% Syntax:
%    checkDynParameter(field,sys,params,options,listname)
%
% Inputs:
%    field - struct field in params / options
%    sys - object of system class
%    func - function
%    params - struct containing model parameters
%    options - struct containing algorithm parameters
%    listname - 'params' or 'options'
%
% Outputs:
%    failedChecks - failed checks
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: checkDynParameterParams, checkDynParameterOptions

% Authors:       Tobias Ladner
% Written:       05-October-2023
% Last update:   09-October-2023 (TL, split options/params)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. init checks
checks = struct('checkfun',{},'errorId',{});

% 2. split list for params / options for more transparency / readability
switch listname
    case 'params'
        % search for checks in params
        checks = checkDynParameterParams(field,sys,func,params,options,checks);

    case 'options'
        % search for checks in options
        checks = checkDynParameterOptions(field,sys,func,params,options,checks);
    
    otherwise
        % unknown listname
        throw(CORAerror('CORA:specialError','listname has to be ''params'' or ''options''.'))

end

% 3. check parameters
failedChecks = aux_checkParams(listname,field,sys,params,options,checks);

end


% Auxiliary functions -----------------------------------------------------

function failedChecks = aux_checkParams(listname,field,sys,params,options,checks)

% init indices for which a check has failed (note: the '[]' is crucial for
% the concatenation later on!)
failedChecks = struct([]);
failedIdx = 0;

% special name-reading for parameters which are part of a group
% of parameters organized in a struct
if strcmp(listname,'params')
    list = params;
elseif  strcmp(listname,'options')
    list = options;
else
    % unknown listname
    throw(CORAerror('CORA:specialError','listname has to be ''params'' or ''options''.'))
end
if contains(field,'.')
    dotIdx = strfind(field,'.');
    firstname = field(1:dotIdx-1);
    secondname = field(dotIdx+1:end);
    val = list.(firstname).(secondname);
else
    val = list.(field);
end

% loop over all check functions for the given parameter
for i=1:length(checks)
    
    % read out identifier for error message
    errmsgid = checks(i).errorId;

    % check function call
    if isempty(errmsgid) 
        % c_* validation functions return also msg
        [rescheck,msg] = checks(i).checkfun(val);
    else 
        % standard case
        rescheck = checks(i).checkfun(val);
        msg = getErrorMessage(errmsgid);
    end

    if ~rescheck
        % check was not ok

        if VALIDATEOPTIONS_ERRORS
            % error message if validation failed
            throw(CORAerror('CORA:specialError', ...
            sprintf('%s.%s: %s.',listname, field, msg)));
        else
            % append to list of failed checks
            failedIdx = failedIdx + 1;
            failedChecks(failedIdx).parameter = blueprint.name{i};
            failedChecks(failedIdx).message = msg;            
        end
    end
end

end

% ------------------------------ END OF CODE ------------------------------
