function res = contains(S1,S2,varargin)
% contains - determines if a set contains another set or a point
%
% Syntax:
%    res = contains(S1,S2)
%    res = contains(S1,S2,method)
%    res = contains(S1,S2,method,tol)
%    res = contains(S1,S2,method,tol,maxEval)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object, numeric array
%    method - method for computation ('exact' or 'approx')
%    tol - tolerance
%    maxEval - see zonotope/contains
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       18-August-2022
% Last update:   23-November-2022 (MW, add classname as input argument)
% Last revision: 27-March-2023 (MW, restructure relation to subclass)

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
if nargin < 2
    throw(CORAerror('CORA:notEnoughInputArgs',2));
elseif nargin > 5
    throw(CORAerror('CORA:tooManyInputArgs',5));
end

% right order of objects
% [S1,S2] = findClassArg(S1,S2,classname);

% parse input arguments (dummy value for maxEval)
[type,tol,maxEval] = setDefaultValues({'exact',100*eps,0},varargin); 

if isa(S1,'zonotope')
    % check input arguments
    inputArgsCheck({{S1,'att','contSet','scalar'};
                    {S2,'att',{'contSet','numeric'}};
                    {type,'str',{'exact','approx','venum','polymax','opt','st'}};
                    {tol,'att','numeric',{'scalar','nonnegative','nonnan'}};
                    {maxEval,'att','numeric',{'scalar','nonnegative','nonnan'}}});

    % default value for maxEval depends on in-body zonotope
    if isa(S2,'zonotope') && strcmp(type,'opt') && maxEval == 0
        m1 = size(generators(S2),2);
        maxEval = max(500,200*m1);
    end

elseif isa(S1,'zonoBundle')
    % more types...
    inputArgsCheck({{S1,'att','contSet'};
                    {S2,'att',{'contSet','numeric'}};
                    {type,'str',{'exact','approx','exact:zonotope','exact:conZonotope','exact:polytope'}};
                    {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});

else
    % check input arguments
    inputArgsCheck({{S1,'att','contSet'};
                    {S2,'att',{'contSet','numeric'}};
                    {type,'str',{'exact','approx'}};
                    {tol,'att','numeric',{'scalar','nonnegative','nonnan'}}});
    
end


% check dimension mismatch (has to come after empty set case!)
equalDimCheck(S1,S2);

% call subclass method
% try
    res = contains_(S1,S2,type,tol,maxEval);

% catch ME
%     % empty set cases
% 
%     % inner-body is point = [] or empty set of any contSet class
%     if representsa_(S2,'emptySet',1e-8)
%         res = true;
% 
%     elseif representsa_(S1,'emptySet',1e-8)
%         % outer body is empty: containment would only be fulfilled if inner
%         % body is empty too, which is handled above
%         res = false;
%         
%     else
%         rethrow(ME);
%     end
% end

% ------------------------------ END OF CODE ------------------------------
