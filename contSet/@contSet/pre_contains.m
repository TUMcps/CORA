function [res,vars] = pre_contains(classname,S1,S2,varargin)
% pre_contains - pre-processing for contSet/contains functions including
%    1. setting of default values
%    2. checking of input arguments
%    3. computation of special cases (e.g., empty set)
%    the class name of the calling function is handed as an input argument
%    as a read-out from the call stack would be computationally inefficient
%
% Syntax:
%    [res,vars] = pre_contains(classname,S1,S2)
%    [res,vars] = pre_contains(classname,S1,S2,method)
%    [res,vars] = pre_contains(classname,S1,S2,method,tol)
%    [res,vars] = pre_contains(classname,S1,S2,method,tol,maxEval)
%
% Inputs:
%    classname - class of calling function
%    S1 - contSet object
%    S2 - contSet object, numeric array
%    method - method for computation ('exact' or 'approx')
%    tol - tolerance
%    maxEval - see zonotope/contains
%
% Outputs:
%    res - true/false whether result has been computed here
%    vars - cell-array whose content depends on value of res:
%           res = true: value of contains
%           res = false: checked input argument arguments incl. default
%                        values for vertices function evaluation
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      18-August-2022
% Last update:  23-November-2022 (MW, add classname as input argument)
% Last revision:---

%------------- BEGIN CODE --------------

% check number of input arguments
if nargin < 3
    throw(CORAerror('CORA:notEnoughInputArgs',3));
elseif nargin > 6
    throw(CORAerror('CORA:tooManyInputArgs',6));
end

% right order of objects
[S1,S2] = findClassArg(S1,S2,classname);

% parse input arguments (dummy value for maxEval)
[type,tol,maxEval] = setDefaultValues({'exact',100*eps,0},varargin); 

if strcmp(classname,'zonotope')
    % check input arguments
    inputArgsCheck({{S1,'att',classname,'scalar'};
                    {S2,'att',{'contSet','numeric'}};
                    {type,'str',{'exact','approx','venum','polymax','opt','st'}};
                    {tol,'att','numeric',{'scalar','nonnegative','nonnan'}};
                    {maxEval,'att','numeric',{'scalar','nonnegative','nonnan'}}});

    % default value for maxEval depends on in-body zonotope
    if isa(S2,'zonotope') && strcmp(type,'opt') && maxEval == 0
        m1 = size(generators(S2),2);
        maxEval = max(500,200*m1);
    end

elseif strcmp(classname,'zonoBundle')
    % more types...
    inputArgsCheck({{S1,'att',classname};
                    {S2,'att',{'contSet','numeric'}};
                    {type,'str',{'exact','approx','exact:zonotope','exact:conZonotope','exact:polytope'}};
                    {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

else
    % check input arguments
    inputArgsCheck({{S1,'att',classname};
                    {S2,'att',{'contSet','numeric'}};
                    {type,'str',{'exact','approx'}};
                    {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});
    
end

% empty set cases
if (isnumeric(S2) && isempty(S2)) || (isa(S2,'contSet') && isemptyobject(S2))
    res = true;
    vars{1} = true;
    return
elseif isemptyobject(S1)
    % case S2 also empty is excluded, therefore always false
    res = true;
    vars{1} = false;
    return
end

% check dimension mismatch (has to come after empty set case!)
equalDimCheck(S1,S2);

% result has to be computed in calling function
res = false;
vars = {S1,S2,type,tol,maxEval};

%------------- END OF CODE --------------
