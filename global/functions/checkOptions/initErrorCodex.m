function initErrorCodex
% initErrorCodex - initializes all possible error messages that can occur
%    in validateOptions; the first argument to the auxiliary function
%    appendErrMsg is the identifier that is used in the fourth input
%    argument to add2params/add2options in every config-file
%
% Syntax:
%    initErrorCodex
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      08-February-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% define as global variables for simpler syntax in config file
global codex;

% init empty struct
codex = struct();
% introduce fields
codex.id = {};
codex.text = {};

% add error messages using aux. func. 'appendErrMsg'
% note: write text such that the first word is the name of the param/option
%       -> example: 'has to be a scalar value'
%                   becomes, e.g.,
%                   'options.timeStep has to be a scalar value.'
%       (the dot will also be added automatically)

% no error message
appendErrMsg('none','should be true by default (report issue)');

% generic error message
appendErrMsg('???','... please check the description given in the manual');

% error message in c_* function
appendErrMsg('','see c_* function');

% type checks
appendErrMsg('isscalar','has to be a scalar value');
appendErrMsg('isnumeric','has to be a numeric value');
appendErrMsg('islogical','has to be a logical (true/false) value');
appendErrMsg('ischar','has to be a char array');
appendErrMsg('iscell','has to be a cell array');
appendErrMsg('isvector','has to be a vector');
% isa(...) checks
appendErrMsg('isafunction_handle','has to be a function handle');
appendErrMsg('isacontSet','has to be an object of class contSet');
appendErrMsg('isazonotope','has to be an object of class zonotope');

% absolute comparisons
appendErrMsg('gtzero','has to be a value greater than 0');
appendErrMsg('gezero','has to be a value greater than or equal to 0');
appendErrMsg('vectorgezero','has to have a value greater than or equal to 0 in all entries');
appendErrMsg('geone','has to be a value greater than or equal to 1');
appendErrMsg('vectorgeone','has to have a value greater than or equal to 1 in all entries');
appendErrMsg('normalized','has to be a value in [0,1]');
appendErrMsg('integer','has to be an integer value');
appendErrMsg('integerorInf','has to be an integer value or Inf');
appendErrMsg('vectororinterval','has to be a vector or an object of class interval');
% just tensorOrder
appendErrMsg('2or3','has to be either 2 or 3');


% not in list of admissible values specified in getMembers
% (second input arg to appendErrMsg just dummy)
appendErrMsg('memberR0','...');
appendErrMsg('memberU','...');
appendErrMsg('memberUsim','...');
appendErrMsg('memberV','...');
appendErrMsg('memberW','...');
appendErrMsg('memberlinAlg','...');
appendErrMsg('memberlinAlg4HA','...');
appendErrMsg('memberalg','...');
appendErrMsg('memberalgInner','...');
appendErrMsg('memberalg4param','...');
appendErrMsg('memberalg4observe','...');
appendErrMsg('memberreductionTechnique','...');
appendErrMsg('memberreductionTechnique4nlsys','...');
appendErrMsg('memberreductionTechniqueUnderApprox','...');
appendErrMsg('memberlagrangeRem.simplify','...');
appendErrMsg('memberlagrangeRem.method','...');
appendErrMsg('memberlagrangeRem.zooMethods','...');
appendErrMsg('memberlagrangeRem.optMethod','...');
appendErrMsg('memberrestructureTechnique','...');
appendErrMsg('membercontractor','...');
appendErrMsg('memberguardIntersect','...');
appendErrMsg('memberenclose','...');

% comparisons to obj
appendErrMsg('eqsysdim','has to be equal to the system dimension');
appendErrMsg('eqinput','has to be equal to the number of inputs');
appendErrMsg('eqoutput','has to be equal to the number of outputs');
appendErrMsg('eqconstr','has to be equal to the number of constraints');
appendErrMsg('eqparam','has to be equal to the number of parameters');

% comparisons to other params/options
appendErrMsg('getStart','has to be greater than params.tStart');
appendErrMsg('intsteps','has to divide the time horizon into an integer number steps');
appendErrMsg('letaylorTerms','has to be less or equal to options.taylorTerms');
appendErrMsg('lecomp','has to be less or equal to the number of components');
appendErrMsg('leloc','has to be less or equal to the number of locations');
appendErrMsg('lelocplus1','has to be less or equal to the number of locations plus one');
appendErrMsg('eqreachSteps','has to be equal to the number of time steps');

% note: global variable 'codex' cleared at the end of validateOptions

end

% Auxiliary Function
function appendErrMsg(newid,newtext)

% enable access to codex
global codex;

codex.id{end+1,1} = newid;
% if id contains member, get list of possible values from getMembers.m
if contains(newid,'member')
    codex.text{end+1,1} = ['has to match one of the following:\n' ...
        strjoin(getMembers(erase(newid,'member')),', ')];
else
    codex.text{end+1,1} = newtext;
end

end

%------------- END OF CODE --------------