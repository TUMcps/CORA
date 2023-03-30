classdef halfspace < contSet
% halfspace - object constructor for halfspaces 
%
% Description:
%    This class represents halfspace objects defined as
%    {x | c'*x <= d}.
%
% Syntax:  
%    obj = halfspace()
%    obj = halfspace(hs)
%    obj = halfspace(c,d)
%
% Inputs:
%    c - normal vector of the halfspace (n-by-1)
%    d - distance to the origin (scalar)
%
% Outputs:
%    obj - generated halfspace object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: example_halfspace.m

% Author:       Matthias Althoff
% Written:      06-June-2011
% Last update:  14-Aug-2019
%               02-May-2020 (add property validation)
%               19-March-2021 (MW, error messages)
%               14-December-2022 (TL, property check in inputArgsCheck)
% Last revision:---

%------------- BEGIN CODE --------------


properties (SetAccess = private, GetAccess = public)
    c = [];
    d = 0;
end
    
methods
    %class constructor
    function obj = halfspace(varargin)
        
        % parse input
        if nargin > 2
            throw(CORAerror('CORA:tooManyInputArgs',2));
        end

        if nargin==1
            hs = varargin{1};
            inputArgsCheck({{hs, 'att', 'halfspace'}})
            obj = hs;
            return
        end

        [c, d] = setDefaultValues({[], 0}, varargin);
        

        if nargin > 0
            inputArgsCheck({ ...
                {c, 'att', 'numeric', {'finite', 'vector'}}; ...
                {d, 'att', 'numeric', {'finite', 'scalar'}}; ...
            })
            c = reshape(c, [], 1); % column vector
        end

        % assign properties
        obj.c = c;
        obj.d = d;
    end
end

methods (Static = true)
    hs = generateRandom(varargin) % generates random halfspace
end


end

%------------- END OF CODE --------------