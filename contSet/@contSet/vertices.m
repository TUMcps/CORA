function res = vertices(S,varargin)
% vertices - computes the vertices of a set
%
% Syntax:
%    res = vertices(S)
%    res = vertices(S,method)
%
% Inputs:
%    S - contSet object
%    method - method for computation of vertices
%
% Outputs:
%    res - array of vertices
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       18-August-2022
% Last update:   23-November-2022 (MW, add classname as input argument)
%                12-July-2023 (TL, corrected dimension of empty vertices)
% Last revision: 27-March-2023 (MW, restructure relation to subclass)

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
if nargin < 1
    throw(CORAerror('CORA:notEnoughInputArgs',1));
elseif nargin > 2
    throw(CORAerror('CORA:tooManyInputArgs',2));
end

% default values and input argument check
if isa(S,'polytope')
    method = setDefaultValues({'lcon2vert'},varargin);
    inputArgsCheck({{S,'att','polytope'}, ...
                    {method,'str',{'cdd','lcon2vert','comb'}}});
else
    method = setDefaultValues({'convHull'},varargin);
    inputArgsCheck({{S,'att','contSet'};
                    {method,'str',{'convHull','iterate','polytope'}}});
end

% call subclass method
try
    res = vertices_(S,method);

catch ME
    % catch empty set case
    if representsa_(S,'emptySet',eps)
        res = [];
    else
        rethrow(ME);
    end
end

if isempty(res)
    % create res with proper dimensions
    res = zeros(dim(S),0);
end
    
end

% ------------------------------ END OF CODE ------------------------------
