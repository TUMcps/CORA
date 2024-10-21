function res = vertices(S,varargin)
% vertices - computes the vertices of a set
%
% Syntax:
%    res = vertices(S)
%    res = vertices(S,method,varargin)
%
% Inputs:
%    S - contSet object
%    method - method for computation of vertices
%    varargin - further parameters
%
% Outputs:
%    res - array of vertices
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/polygon

% Authors:       Mark Wetzlinger
% Written:       18-August-2022
% Last update:   23-November-2022 (MW, add classname as input argument)
%                12-July-2023 (TL, corrected dimension of empty vertices)
%                12-July-2024 (MW, disable 'comb' for polytopes)
% Last revision: 27-March-2023 (MW, restructure relation to subclass)

% ------------------------------ BEGIN CODE -------------------------------

% check number of input arguments
narginchk(1,3);

% default values and input argument check
[S,method,addargs] = aux_parseInput(S,varargin{:});

% call subclass method
try
    res = vertices_(S,method,addargs{:});

catch ME
    % catch empty set case
    if representsa_(S,'emptySet',eps,'linearize',0,1)
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


% Auxiliary functions -----------------------------------------------------

function [S,method,addargs] = aux_parseInput(S,varargin)
    if isa(S,'polytope')
        % uses different methods
        method = setDefaultValues({'lcon2vert'},varargin);
        inputArgsCheck({{S,'att','polytope'}, ...
                        {method,'str',{'cdd','lcon2vert'}}});
    elseif isa(S,'conPolyZono')
        % 'method' is number of splits
        method = setDefaultValues({10},varargin);
        inputArgsCheck({{S,'att','conPolyZono'}, ...
                        {method,'att','numeric',{'scalar','nonnan'}}});
    elseif isa(S,'conZonotope')
        [method,numDirs] = setDefaultValues({'default',1},varargin);
        inputArgsCheck({ ...
            {S,'att','conZonotope'}; ...
            {method,'str',{'default','template'}}; ...
            {numDirs,'att','numeric','isscalar'}; ...
        })
        varargin = {method,numDirs};
    
    else
        method = setDefaultValues({'convHull'},varargin);
        inputArgsCheck({{S,'att','contSet'};
                        {method,'str',{'convHull','iterate','polytope'}}});
    end
    addargs = varargin(2:end);
end

% ------------------------------ END OF CODE ------------------------------
