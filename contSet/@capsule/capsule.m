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
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    c (:,1) {mustBeNumeric,mustBeFinite} = []; % center
    g (:,1) {mustBeNumeric,mustBeFinite} = []; % generator
    r (1,1) {mustBeNonnegative,mustBeFinite} = 0; % radius
end
   
methods

    function Obj = capsule(varargin)
        
        % default constructor
        if nargin == 0
            
            Obj.dimension = 0;
            
        else
            
            % If 1 argument is passed
            if nargin == 1
                if isa(varargin{1},'capsule')
                    % copy constructor
                    Obj = varargin{1};
                else
                    % set center
                    Obj.c = varargin{1}; 
                    Obj.g = zeros(length(Obj.c),1);
                end
            elseif nargin == 2
                Obj.c = varargin{1};
                if isscalar(varargin{2})
                    Obj.r = varargin{2};
                    Obj.g = zeros(length(Obj.c),1);
                elseif isvector(varargin{2}) && length(Obj.c) ~= length(varargin{2})
                    % dimension mismatch
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Dimension mismatch between center and generator.'));
                else
                    Obj.g = varargin{2};
                end
            elseif nargin == 3
                % set all values
                Obj.c = varargin{1};
                if isvector(varargin{2}) && length(Obj.c) ~= length(varargin{2})
                    % dimension mismatch
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Dimension mismatch between center and generator.'));
                end
                Obj.g = varargin{2};
                Obj.r = varargin{3};
                
            elseif nargin > 3
                % too many input arguments
                throw(CORAerror('CORA:tooManyInputArgs',3));
            end

            % set parent object properties
            Obj.dimension = length(Obj.c);
        end
    end
end

methods (Static = true)
    C = generateRandom(varargin) % generates random capsule
end

end
%------------- END OF CODE --------------