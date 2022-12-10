function [res,vars] = pre_randPoint(classname,S,varargin)
% pre_randPoint - pre-processing for contSet/randPoint functions including
%    1. setting of default values
%    2. checking of input arguments
%    3. computation of special cases (e.g., empty set)
%    the class name of the calling function is handed as an input argument
%    as a read-out from the call stack would be computationally inefficient
%
% Syntax:
%    [res,vars] = pre_randPoint(classname,S)
%    [res,vars] = pre_randPoint(classname,S,N)
%    [res,vars] = pre_randPoint(classname,S,N,type)
%    [res,vars] = pre_randPoint(classname,S,N,type,pr)
%
% Inputs:
%    classname - class of calling function
%    S - contSet object
%    N - (optional) number of samples (integer) or 'all' (for sets with
%                   vertices)
%    type - (optional) sampling method ('standard', 'extreme', 'gaussian')
%    pr - (optional) probability that sampled point is within the ellipsoid
%                    outer-approximation of S (only type = 'gaussian')
%
% Outputs:
%    res - true/false whether result has been computed here
%    vars - cell-array whose content depends on value of res:
%           res = true: value of randPoint
%           res = false: checked input argument arguments incl. default
%                        values for randPoint function evaluation
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      19-August-2022
% Last update:  23-November-2022 (MW, add classname as input argument)
% Last revision:---

%------------- BEGIN CODE --------------

% check number of input arguments
if nargin < 1
    throw(CORAerror('CORA:notEnoughInputArgs',1));
elseif nargin > 5
    throw(CORAerror('CORA:tooManyInputArgs',5));
end

% set default values for number of samples, method, and probability
[N,type,pr] = setDefaultValues({1,'standard',0.7},varargin{:});

% N can be numeric or 'all'
if isnumeric(N)
    checkN = {N,'att','numeric',{'scalar','integer','positive'}};
else
    checkN = {N,'str','all'};
end
% check input arguments
inputArgsCheck({{S,'att',classname};
                 checkN; % see above...
                {type,'str',{'standard','extreme','gaussian'}};
                {pr,'att','numeric',{'<=',1,'>=',0}}});

% if N = 'all', then type has to be 'extreme'
if ischar(N) && strcmp(N,'all') && ~strcmp(type,'extreme')
    throw(CORAerror('CORA:wrongValue','third',...
        "If the number of points is 'all', the type has to be 'extreme'."));
end

% special case: contSet is empty set
if isemptyobject(S)
    res = true;
    vars{1} = []; return
end

% result has to be computed in calling function
res = false;
vars = {S,N,type,pr};

%------------- END OF CODE --------------
