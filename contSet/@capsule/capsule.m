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
%    obj = capsule(c,r)
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
%     
%    C = capsule(c,g,r);
%
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
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    c; % center
    g; % generator
    r; % radius
end
   
methods

    function obj = capsule(varargin)
        if nargin > 3
            % too many input arguments
            throw(CORAerror('CORA:tooManyInputArgs',3));
        end

        % parse input
        if nargin == 1 && isa(varargin{1},'capsule')
            % copy constructor
             obj = varargin{1};
             return;
        end
        
        if nargin == 0
            obj.c = [];
            obj.g = [];
            obj.r = 0;
            
        else
            c = varargin{1};
            g = zeros(length(c),1);
            r = 0;

            if nargin == 2
                v2 = varargin{2};
                if isscalar(v2)
                    r = v2;
                elseif isvector(v2) && length(c) == length(v2)
                    g = v2;
                elseif ~isempty(v2)
                    % dimension mismatch
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Dimension mismatch between center and generator.'));
                end
            end

            if nargin == 3
               g = varargin{2};
               r = varargin{3};
            end

            inputArgsCheck({ ...
                {c, 'att', {'numeric','numeric'}, {{'nonempty','finite','column'},{'size',[]}} }; ...
                {g, 'att', {'numeric','numeric'}, {{'nonempty','finite','column'},{'size',[]}} }; ...
                {r, 'att', 'numeric', {'finite', 'nonnegative', 'scalar'}};
            })

            % set properties
            obj.c = c;
            obj.g = g;
            obj.r = r;
        end
    end
end

methods (Static = true)
    C = generateRandom(varargin) % generates random capsule
end

end
%------------- END OF CODE --------------