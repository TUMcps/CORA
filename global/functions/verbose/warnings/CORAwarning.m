function CORAwarning(identifier, varargin)
% CORAwarning - prints a CORA warning into the command line.
%    Can be controlled via the CORA_WARNINGS_ENABLED macro
%
% Syntax:
%    CORAwarning(identifier, varargin)
%
% Inputs:
%    identifier - char
%    varargin - depends on identifier, see below
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: CORA_WARNINGS_ENABLED, CORAerror

% Authors:       Tobias Ladner
% Written:       17-April-2024
% Last update:   07-October-2024 (MW, add CORA:interface)
%                17-September-2025 (TL, made call stack less verbose)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(2,Inf);
aux_checkIdentifier(identifier)

% check if CORA warning should be shown (global/specific)
warnState = warning();
if (strcmp(warnState(1).identifier,'all') && strcmp(warnState(1).state,'off')) ...
    || ~CORA_WARNINGS_ENABLED || ~CORA_WARNINGS_ENABLED(identifier)
    % do not show warning
    return
end

% build specific warning text
switch identifier
    case 'CORA:deprecated'
        % identifier, type, name, version, replacement, reason

        % read input
        narginchk(6,6);
        type = varargin{1};
        name = varargin{2};
        version = varargin{3};
        replacement = varargin{4};
        reason = varargin{5};

        % compose replacement and reason
        replacement = aux_formatString(replacement);
        reason = aux_formatString(reason);
        
        desc = sprintf('The %s ''%s'' is deprecated (since %s) and will be removed in a future release.\n  %s\n  %s', ...
            type, name, version, replacement, reason);

    case 'CORA:interface'
        % read input
        narginchk(3,3);
        name = varargin{1};
        version = varargin{2};
        
        desc = sprintf('The interface of the function ''%s'' has changed (since %s):\n  Please check the function header for details', ...
            name, version);

    otherwise
        % identifier, msg, args for sprintf
        desc = varargin{1};
        desc = aux_formatString(desc);
        desc = sprintf(desc,varargin{2:end});
end

% show warning (do not use warning to avoid showing full stack trace all the time)

% read out stacktrace
stacktrace = evalc('dbstack');
% remove CORAwarning line and split remaining lines
stacktrace = split(stacktrace,newline);
stacktrace = stacktrace(2:end);
if numel(stacktrace) == 1 && isempty(stacktrace{1})
    % directly running in command window does not have a stack trace
    stacktrace{1} = 'Run from Command Window';
end
% join and escape full stack trace
fullstack = strjoin(stacktrace,'\n ');
fullstack = escapewarning(fullstack);

% show warning
fprintf('[\b<strong>CORA warning:</strong> %s\n %s (<a href="matlab:fprintf(descapewarning(''Full call stack:\\n %s\\n''))">show full call stack</a>).]\b\n', ...
    desc, stacktrace{1}, fullstack)

end


% Auxiliary functions -----------------------------------------------------

function aux_checkIdentifier(identifier)
    % test identifier 
    % (please also add new identifiers to test_cora_warnings_enabled)
    inputArgsCheck({{identifier,'str',{ ...
        'CORA:app', ...
        'CORA:contDynamics', ...
        'CORA:contSet', ...
        'CORA:converter', ...
        'CORA:discDynamics', ...
        'CORA:examples', ...
        'CORA:global', ...
        'CORA:hybridDynamics', ...
        'CORA:manual', ...
        'CORA:matrixSets', ...
        'CORA:models', ...
        'CORA:nn', ...
        'CORA:specification', ...
        'CORA:unitTests', ...
         ... % special warnings
        'CORA:plot', ...
        'CORA:solver', ...
        'CORA:deprecated', ...
        'CORA:redundant', ...
        'CORA:interface', ...
    }}})
end

function str = aux_formatString(str)
    % format line breaks
    str = replace(str,'\n','\n  ');
    str = compose(str);
    str = str{1};
end

% ------------------------------ END OF CODE ------------------------------
