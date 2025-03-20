classdef contSet
% contSet - abstract superclass for continuous sets
%
% Syntax:
%    S = contSet()
%    S = contSet(S)
%
% Inputs:
%    S - contSet object
%
% Outputs:
%    obj - generated contSet object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       02-May-2007 
% Last update:   04-May-2020 (MW, transition to classdef)
%                01-June-2022 (MW, add CORAerror)
%                22-March-2023 (MW, remove deprecated property dimension)
%                22-September-2024 (MW, add hidden property)
% Last revision: 16-June-2023 (MW, restructure using standardized workflow)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = protected, GetAccess = public)
    % no properties
end

properties (SetAccess = protected, Hidden)
    % ordering of set representation (roughly along their generality)
    % -> lower integer value implements commutative set operation
    % 0 emptySet
    % 10 fullspace
    % 15 polygon
    % 20 levelSet
    % 30 conPolyZono
    % 40 spectraShadow
    % 50 ellipsoid
    % 60 capsule
    % 70 polyZonotope
    % 80 polytope
    % 90 conZonotope
    % 100 zonoBundle
    % 110 zonotope
    % 120 interval
    % 1000 contSet
    % 'affine', 'taylm', 'probZonotope', 'zoo' skipped since they do not
    % mingle with the other set representations
    precedence;
    % affected operations: and, convHull, isequal, isIntersecting, or, plus
end

methods

    function obj = contSet(varargin)

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'contSet')
            obj = varargin{1}; return
        end
        assertNarginConstructor(0,nargin);

        % 2. assign properties
        obj.precedence = 1000;
    end

    % prevent instantiation of class arrays
    % -> deleted horizontal concatenation (horzcat), vertical concatenation
    %    (vertcat), concatenation along different dimensions (cat),
    %    assignment to index with () or {} (subsasgn)
    % note: subsasgn is also used for dot access, e.g., Z.c
    % note: read from index (subsref) is not harmful
    % vertical concatenation
    function obj = vertcat(varargin)
        throw(CORAerror('CORA:notSupported',...
            'Given subclass of contSet does not support class arrays.'));
    end
    % horizontal concatenation
    function obj = horzcat(varargin)
        throw(CORAerror('CORA:notSupported',...
            'Given subclass of contSet does not support class arrays.'));
    end
    % concatenation along different dimensions
    function obj = cat(varargin)
        throw(CORAerror('CORA:notSupported',...
            'Given subclass of contSet does not support class arrays.'));
    end
    % assignment to index
    function obj = subsasgn(S,strct,val)
        if ~strcmp(strct(1).type,'.')
            throw(CORAerror('CORA:notSupported',...
                'Given subclass of contSet does not support class arrays.'));
        end
        obj = builtin('subsasgn',S,strct,val);
    end
end

% Static methods
methods (Static = true)
    S = enclosePoints(varargin) % encloses a point cloud with a set
    S = generateRandom(varargin) % generates a random contSet
    S = initEmptySet(type) % instantiates an empty set of a contSet class
    S = empty(type) % instantiates an empty set of a contSet class
    S = Inf(type) % instantiates a fullspace set of a contSet class
end

% Protected access
methods (Access = protected)
    [empty,res,S_conv] = representsa_emptyObject(S,type)
    [abbrev,printOrder] = getPrintSetInfo(S)
end

% plotting
methods (Access = {?contSet})
    han = plot1D(S,varargin);
    han = plot2D(S,varargin);
    han = plot3D(S,varargin);
end

end

% ------------------------------ END OF CODE ------------------------------
