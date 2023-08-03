function S = generateRandom(varargin)
% generateRandom - Generates a random contSet
%
% Syntax:
%    S = contSet.generateRandom()
%    S = contSet.generateRandom('Dimension',n)
%    S = contSet.generateRandom({@interval, @zonotope}, ...)
%
% Inputs:
%    admissibleSets - cell array of addmissible sets
%    Name-Value pairs (all options, arbitrary order):
%       <'Dimension',n> - dimension
%       ... - further pairs for subclasses
%
% Outputs:
%    S - contSet
%
% Example:
%    S = contSet.generateRandom('Dimension',2);
%    plot(S);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Tobias Ladner
% Written:      05-April-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% parse input
if nargin < 1 || ~iscell(varargin{1})
    admissibleSets = {
        @capsule; ...
        @conPolyZono; ...
        @conHyperplane; ...
        @conZonotope; ...
        @ellipsoid; ...
        @emptySet; ...
        @fullspace; ...
        @halfspace; ...
        @interval; ...
        @levelSet; ...
        @mptPolytope; ...
        @polyZonotope; ...
        @probZonotope; ...
        @zonoBundle; ...
        @zonotope; ...
        };
else
    admissibleSets = varargin{1};
    varargin = varargin(2:end);
end

if CHECKS_ENABLED
    inputArgsCheck({{admissibleSets, 'att', 'cell'}})
    if ~all(cellfun(@(setHandle) isa(setHandle, 'function_handle'), admissibleSets))
        throw(CORAerror("CORA:wrongValue", "first", ...
            ['Admissible sets should be a cell array of contSet function handles, ', ...
            'for which generateRandom can be called.']))
    end
end

set = admissibleSets{randi(length(admissibleSets))};
S = set().generateRandom(varargin{:});

%------------- END OF CODE --------------