function res = isIntersecting(S1,S2,varargin)
% isIntersecting - checks if two sets intersect
%
% Syntax:
%    res = isIntersecting(S1,S2)
%    res = isIntersecting(S1,S2,type)
%
% Inputs:
%    S1,S2 - contSet object
%    type - (optional) type ('exact','approx')
%
% Outputs:
%    res - true/false whether result has been computed here
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       17-August-2022
% Last update:   23-November-2022 (MW, add classname as input argument)
% Last revision: 27-March-2023 (MW, restructure relation to subclass)

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
if nargin < 2
    throw(CORAerror('CORA:notEnoughInputArgs',2));
elseif nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% set default values: type ('exact')
type = setDefaultValues({'exact'},varargin);

% right order of objects
% [S1,S2] = findClassArg(S1,S2,classname);

% check input arguments: two tries depending on value for type
inputArgsCheck({{S1,'att','contSet'},...
                {S2,'att',{'contSet','numeric'}},...
                {type,'str',{'exact','approx'}}});

% check dimension mismatch
equalDimCheck(S1,S2);

% call subclass method
try
    res = isIntersecting_(S1,S2,type);
catch ME
    % empty set case
    if representsa_(S1,'emptySet',1e-8) || representsa_(S2,'emptySet',1e-8)
        res = false;
    else
        rethrow(ME);
    end
end

% ------------------------------ END OF CODE ------------------------------
