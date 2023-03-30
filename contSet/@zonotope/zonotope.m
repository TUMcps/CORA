classdef (InferiorClasses = {?intervalMatrix, ?matZonotope}) zonotope < contSet
% zonotope - object constructor for zonotope objects
%
% Description:
%    This class represents zonotopes objects defined as
%    {c + \sum_{i=1}^p beta_i * g^(i) | beta_i \in [-1,1]}.
%
% Syntax:
%    obj = zonotope()
%    obj = zonotope(c,G)
%    obj = zonotope(Z)
%
% Inputs:
%    c - center vector
%    G - generator matrix
%    Z - center vector and generator matrix Z = [c,G]
%
% Outputs:
%    obj - generated zonotope object
%
% Example: 
%    c = [1;1];
%    G = [1 1 1; 1 -1 0];
%    Z = zonotope(c,G);
%    plot(Z);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      14-September-2006 
% Last update:  22-March-2007
%               04-June-2010
%               08-February-2011
%               18-November-2015
%               05-December-2017 (DG) class is redefined in compliance with
%               the new standard.
%               28-April-2019 code shortened
%               1-May-2020 (NK) new constructor + removed orientation prop.
%               14-December-2022 (TL, property check in inputArgsCheck)
%               29-March-2023 (TL: optimized constructor)
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    % zonotope center and generator Z = [c,g_1,...,g_p]
    Z;

    % halfspace representation of the zonotope
    halfspace;
end

methods

    function obj = zonotope(varargin)

        % parse input
        switch nargin
            case 0
                Z = [];
            case 1
                if isa(varargin{1}, 'zonotope')
                    % copy constructor
                    obj = varargin{1};
                    return;
                end
                Z = varargin{1};
            case 2
                c = varargin{1};
                G = varargin{2};
    
                if CHECKS_ENABLED
                    inputArgsCheck({ ...
                        {c, 'att', 'numeric', 'nonnan'}; ...
                        {G, 'att', 'numeric', 'nonnan'}; ...
                    })
        
                    % check dimensions
                    if isempty(c) && ~isempty(G)
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'Center is empty.'));
                    elseif ~isvector(c)
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'Center is not a vector.'));
                    elseif ~isempty(G) && size(c,1) ~= size(G,1)
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'Dimension mismatch between center and generator matrix.'));  
                    end
                end
                Z = [c, G];
            otherwise
                throw(CORAerror("CORA:tooManyInputArgs", 2))
        end

        inputArgsCheck({{Z, 'att', 'numeric', 'nonnan'}})

        % assign properties
        obj.Z = Z;
        obj.halfspace = [];
    end
end

methods (Static = true)
    Z = generateRandom(varargin) % generate random zonotope
    Z = enclosePoints(points,varargin) % enclose point cloud with zonotope
end

end

%------------- END OF CODE --------------