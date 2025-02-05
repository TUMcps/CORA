function [res,cert,scaling] = contains(S1,S2,varargin)
% contains - determines if a set contains another set or a point
%
% Syntax:
%    res = contains(S1,S2)
%    res = contains(S1,S2,method)
%    res = contains(S1,S2,method,tol)
%    res = contains(S1,S2,method,tol,maxEval)
%    [res,cert] = contains(___)
%    [res,cert,scaling] = contains(___)
%
% Inputs:
%    S1 - contSet object
%    S2 - contSet object, numeric array
%    method - method for computation ('exact' or 'approx', or any method
%           name specific to a certain set representation; see the
%           documentation of the corresponding contains_ function, e.g.,
%           zonotope/contains_ )
%    tol - tolerance
%    maxEval - maximal number of iterations for optimization-based methods.
%            See the corresponding contains_ (e.g., zonotope/contains_)
%            function for more details
%
% Outputs:
%    res - true if S2 is contained in S1, false if S2 is not contained
%          S1, or if the containment could not be certified
%    cert - (optional) returns true iff the result of res could be
%            verified. For example, if res=false and cert=true, S2 is
%            guaranteed to not be contained in S1, whereas if res=false and
%            cert=false, nothing can be deduced (S2 could still be
%            contained in S1).
%            If res=true, then cert=true.
%            Note that computing this certification may marginally increase
%            the runtime.
%    scaling - (optional) assuming S1 has a center, returns the smallest
%            number 'scaling', such that
%            scaling*(S1 - S1.center) + S1.center contains S2.
%            For methods other than 'exact', this may be a lower or upper
%            bound on the exact number. See the corresponding contains_
%            (e.g., zonotope/contains_) for more details.
%            Note that computing this scaling factor may significantly
%            increase the runtime.
%    
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
narginchk(2,5);

% right order of objects
% [S1,S2] = findClassArg(S1,S2,classname);

% parse input arguments (dummy value for maxEval)
[method,tol,maxEval] = setDefaultValues({'exact',100*eps,0},varargin);
% default value for maxEval depends on in-body zonotope
if isa(S2, 'zonotope') && maxEval == 0
    maxEval = max(500, 200*size(generators(S2),2));
elseif maxEval == 0
    maxEval = 200;
end

% check input arguments
inputArgsCheck({{S1,'att','contSet'};
    {S2,'att',{'contSet','taylm','numeric'}};
    {method,'str',{'exact', 'exact:venum', 'exact:polymax',... % Typical exact containment checks
                        'exact:zonotope', 'exact:polytope',... % EXCLUSIVELY for zonotope bundles
                   'approx', 'approx:st', 'approx:stDual',... % Approx methods
                   'opt', ... % Optimization-based methods
                   'sampling', 'sampling:primal', 'sampling:dual',... % Sampling-based methods
                    }};
    {tol,'att','numeric',{'scalar','nonnegative','nonnan'}};
    {maxEval,'att','numeric',{'scalar','nonnegative','nonnan'}}});

% check if the user wants cert and scaling
certToggle = false; scalingToggle = false;
if nargout >= 2
    certToggle = true;
    if nargout >= 3
        scalingToggle = true;
    end
end

% call subclass method
try
    % check dimension mismatch
    equalDimCheck(S1,S2);
    [res,cert,scaling] = contains_(S1,S2,method,tol,maxEval,certToggle,scalingToggle);

catch ME
    % empty set cases

    % inner-body is point = [] or empty set of any contSet class
    if (isnumeric(S2) && isempty(S2)) || (isa(S2,'contSet') && isemptyobject(S2))
        res = true;
        cert = true;
        scaling = 0;

    elseif isemptyobject(S1)
        % outer body is empty: containment would only be fulfilled if inner
        % body is empty too, which is handled above
        res = false;
        cert = true;
        scaling = Inf;
        
    else
        rethrow(ME);
    end
end


% ------------------------------ END OF CODE ------------------------------
