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

% Authors:       Tobias Ladner
% Written:       05-April-2023
% Last update:   09-January-2024 (changed admissibleSets to cell arrays of strings)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
[admissibleSets,varargin] = aux_parseInput(varargin);

% randomly select set
set = admissibleSets{randi(length(admissibleSets))};

% generate random set
S = feval(sprintf('%s.generateRandom',set),varargin{:});

end


% Auxiliary functions -----------------------------------------------------

function [admissibleSets,args] = aux_parseInput(args)
    
    % parse input
    if isempty(args) || ~iscell(args{1})
        admissibleSets = {
            'capsule'; ...
            'conPolyZono'; ...
            % 'conHyperplane'; ...
            'conZonotope'; ...
            'ellipsoid'; ...
            'emptySet'; ...
            'fullspace'; ...
            'halfspace'; ...
            'interval'; ...
            'levelSet'; ...
            'polytope'; ...
            'polyZonotope'; ...
            'probZonotope'; ...
            'zonoBundle'; ...
            'zonotope'; ...
            };
    else
        admissibleSets = args{1};
        args = args(2:end);
    end
    
    % parse input
    if CHECKS_ENABLED
        inputArgsCheck({{admissibleSets, 'att', 'cell'}})
        if ~all(cellfun(@(set) ischar(set) || isstring(set), admissibleSets))
            throw(CORAerror("CORA:wrongValue", "first", ...
                ['Admissible sets should be a cell array of contSet strings, ', ...
                'for which generateRandom can be called.']))
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
