function [res,vars] = pre_isIntersecting(classname,S1,S2,varargin)
% pre_isIntersecting - pre-processing for contSet/isIntersecting functions
%    including
%    1. setting of default values
%    2. checking of input arguments
%    3. computation of special cases (e.g., empty set)
%    the class name of the calling function is handed as an input argument
%    as a read-out from the call stack would be computationally inefficient
%
% Syntax:
%    [res,vars] = pre_isIntersecting(classname,S1,S2)
%    [res,vars] = pre_isIntersecting(classname,S1,S2,type)
%
% Inputs:
%    classname - class of calling function
%    S1,S2 - contSet object
%    type - (optional) type ('exact','approx')
%
% Outputs:
%    res - true/false whether result has been computed here
%    vars - cell-array whose content depends on value of res:
%           res = true: value of isIntersecting
%           res = false: checked input argument arguments incl. default
%                        values for isIntersecting function evaluation
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      17-August-2022
% Last update:  23-November-2022 (MW, add classname as input argument)
% Last revision:---

%------------- BEGIN CODE --------------

% check number of input arguments
if nargin < 2
    throw(CORAerror('CORA:notEnoughInputArgs',2));
elseif nargin > 4
    throw(CORAerror('CORA:tooManyInputArgs',4));
end

% set default values: type ('exact')
type = setDefaultValues({'exact'},varargin);

% right order of objects
[S1,S2] = findClassArg(S1,S2,classname);

% check input arguments: two tries depending on value for type
inputArgsCheck({{S1,'att',classname},...
                {S2,'att','contSet'},...
                {type,'str',{'exact','approx'}}});

% initialize that no result has been computed
res = false;

% special case: contSet is empty set
if isemptyobject(S1) || isemptyobject(S2)
    res = true;
    vars{1} = false; return
end

% check dimension mismatch (has to come after empty set case!)
equalDimCheck(S1,S2);

% result has to be computed in calling function
vars{1} = S1;
vars{2} = S2;
vars{3} = type;

%------------- END OF CODE --------------
