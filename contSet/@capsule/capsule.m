classdef capsule < contSet
% capsule - object constructor for capsules
%
% Description:
%    This class represents capsule objects defined as
%    C := L + S, L = {c + g*a | a \in [-1,1]}, S = {x | ||x||_2 <= r}.
%
% Syntax:  
%    obj = capsule(c)
%    obj = capsule(c,g)
%    obj = capsule(c,g,r)
%
% Inputs:
%    c - center
%    g - generator
%    r - radius
%
% Outputs:
%    obj - capsule object
%
% Example:
%    c = [1;2];
%    g = [2;1];
%    r = 1;
%    C = capsule(c,g,r);
%    plot(C);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Author:       Matthias Althoff
% Written:      04-March-2019
% Last update:  02-May-2020 (MW, add property validation)
%               19-March-2021 (MW, error messages, remove capsule(r) case)
%               14-December-2022 (TL, property check in inputArgsCheck)
% Last revision:16-June-2023 (MW, restructure using auxiliary functions)

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    c; % center
    g; % generator
    r; % radius
end
   
methods

    function obj = capsule(varargin)

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'capsule')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [c,g,r] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(c,g,r,nargin);

        % 4. assign properties
        obj.c = c;
        obj.g = g;
        obj.r = r;

    end
end

methods (Static = true)
    C = enclosePoints(varargin) % enclose point cloud by capsule
    C = generateRandom(varargin) % generates random capsule
end

end


% Auxiliary functions -----------------------------------------------------
function [c,g,r] = aux_parseInputArgs(varargin)
% parse input arguments from user and assign to variables

    % check number of input arguments
    if nargin > 3
        throw(CORAerror('CORA:tooManyInputArgs',3));
    end

    % no input arguments
    if nargin == 0
        c = []; g = []; r = 0;
        return
    end

    % assign center
    c = varargin{1};
    % default values for generator and radius
    g = zeros(length(c),1);
    r = 0;

    % assign generator if given
    if nargin > 1 && ~isempty(varargin{2})
        g = varargin{2};
    end
    % assign radius if given
    if nargin > 2
        r = varargin{3};
    end

end

function aux_checkInputArgs(c,g,r,n_in)
% check correctness of input arguments

    % only check if macro set to true
    if CHECKS_ENABLED && n_in > 0

        inputArgsCheck({ ...
            {c, 'att', {'numeric','numeric'}, {{'nonempty','finite','column'},{'size',[]}} }; ...
            {g, 'att', {'numeric','numeric'}, {{'nonempty','finite','column'},{'size',[]}} }; ...
            {r, 'att', 'numeric', {'finite', 'nonnegative', 'scalar'}};
        })

        % dimension of center and generator must match
        if ~all(size(c) == size(g))
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Dimension of center and generator do not match.'));
        end

    end

end

%------------- END OF CODE --------------