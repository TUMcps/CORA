function msg = getErrorMessage(id)
% getErrorMessage - returns the error message corresponding to id
%
% Syntax:
%    msg = getErrorMessage(id)
%
% Inputs:
%    id - identifier for error message
%
% Outputs:
%    msg - error message
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       26-May-2021
% Last update:   19-May-2023 (MW, fix for unit test)
% Last revision: 19-June-2023 (MW, combine with initErrorCodex)

% ------------------------------ BEGIN CODE -------------------------------

% input preprocessing
inputArgsCheck({{id,'att','char'}});

switch id
    
    case 'none'
        % no error message
        msg = 'should be true by default (report issue)';

    case '???'
        % generic error message
        msg = '... please check the description given in the manual';

    case ''
        % error message in c_* function
        msg = 'see c_* function';

    % type checks
    case 'isscalar'
        msg = 'has to be a scalar value';
    case 'isnumeric'
        msg = 'has to be a numeric value';
    case 'islogical'
        msg = 'has to be a logical (true/false) value';
    case 'ischar'
        msg = 'has to be a char array';
    case 'iscell'
        msg = 'has to be a cell array';
    case 'isvector'
        msg = 'has to be a vector';
    case 'isstruct'
        msg = 'has to be a struct';
    case 'isafunction_handle'
        msg = 'has to be a function handle';
    case 'isacontSet'
        msg = 'has to be an object of class contSet';
    case 'isazonotope'
        msg = 'has to be an object of class zonotope';
    case 'isareachSet'
        msg = 'has to be an object of class reachSet';
    case 'istestCase'
        msg = 'has to be an object of class testCase';

    % absolute comparisons
    case 'gtzero'
        msg = 'has to be a value greater than 0';
    case 'gezero'
        msg = 'has to be a value greater than or equal to 0';
    case 'vectorgezero'
        msg = 'has to have a value greater than or equal to 0 in all entries';
    case 'geone'
        msg = 'has to be a value greater than or equal to 1';
    case 'vectorgeone'
        msg = 'has to have a value greater than or equal to 1 in all entries';
    case 'normalized'
        msg = 'has to be a value in [0,1]';
    case 'integer'
        msg ='has to be an integer value';
    case 'integerorInf'
        msg = 'has to be an integer value or Inf';
    case 'vectororinterval'
        msg ='has to be a vector or an object of class interval';
    case 'startatzero'
        msg = 'has to start at 0';
    case '2or3'
        % just tensorOrder
        msg = 'has to be either 2 or 3';

    % not in list of admissible values specified in getMembers
    % (second input arg to appendErrMsg just dummy)
    case {'memberR0','memberR0conf','memberU','memberUsim','memberV','memberW',...
            'membersafeSet','memberunsafeSet','memberlinAlg',...
            'memberlinAlg4HA','memberalg','memberalgInner',...
            'memberalg4DT','memberalg4DA','memberalg4param',...
            'memberalg4observe','memberreductionTechnique',...
            'memberreductionTechnique4nlsys',...
            'memberreductionTechniqueUnderApprox',...
            'memberlagrangeRem.simplify','memberlagrangeRem.method',...
            'memberlagrangeRem.zooMethods','memberlagrangeRem.optMethod',...
            'memberrestructureTechnique','membercontractor',...
            'memberguardIntersect','memberenclose','membertype',...
            'membernorm','memberconfAlgSynth','memberconfAlgCheck',...
            'memberreachAlg','memberarmaxAlg'}
        msg = ['has to match one of the following:\n      ' ...
            strjoin(getMembers(erase(id,'member')),', ')];

    % comparisons to obj
    case 'eqsysdim'
        msg = 'has to be equal to the system dimension';
    case 'eqinput'
        msg ='has to be equal to the number of inputs';
    case 'eqoutput'
        msg ='has to be equal to the number of outputs';
    case 'eqconstr'
        msg = 'has to be equal to the number of constraints';
    case 'eqparam'
        msg = 'has to be equal to the number of parameters';

    % comparisons to other params/options
    case 'getStart'
        msg = 'has to be greater than params.tStart';
    case 'idx1eqtStart'
        msg = 'has to be equal to params.tStart in the first index';
    case 'intsteps'
        msg = 'has to divide the time horizon into an integer number steps';
    case 'letaylorTerms'
        msg = 'has to be less or equal to options.taylorTerms';
    case 'lecomp'
        msg = 'has to be less or equal to the number of components';
    case 'leloc'
        msg = 'has to be less or equal to the number of locations';
    case 'lelocplus1'
        msg = 'has to be less or equal to the number of locations plus one';
    case 'eqreachSteps'
        msg = 'has to be equal to the number of time steps';
    case 'equ'
        msg = 'has to be equal to the length of the time-varying input vector';

    otherwise
        % unknown id
        msg = '';

end

end

% ------------------------------ END OF CODE ------------------------------
