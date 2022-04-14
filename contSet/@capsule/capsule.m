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
%    input1 - center
%    input2 - generator
%    input3 - radius
%
% Outputs:
%    obj - generated capsule object
%
% Example:
%    c = [1;2];
%    g = [2;1];
%    r = 1;
%     
%    C = capsule(c,g,r);
%
%    plot(C);
%    axis equal
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval, polytope

% Author:       Matthias Althoff
% Written:      04-March-2019
% Last update:  02-May-2020 (MW, add property validation)
%               19-March-2021 (MW, errConstructor, remove capsule(r) case)
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
                end
            elseif nargin == 2
                if isvector(varargin{1}) && isvector(varargin{2}) && ...
                    length(varargin{1}) == length(varargin{2})
                    % set center and generator
                    Obj.c = varargin{1}; 
                    Obj.g = varargin{2};
                elseif isvector(varargin{1}) && isscalar(varargin{2})
                    % set center and radius
                    Obj.c = varargin{1}; 
                    Obj.r = varargin{2};
                else
                    % dimension mismatch
                    [id,msg] = errConstructor(); error(id,msg);
                end
                    
            elseif nargin == 3
                % set all values
                Obj.c = varargin{1}; 
                Obj.g = varargin{2};
                Obj.r = varargin{3};
                
            elseif nargin > 3
                % too many input arguments
                [id,msg] = errConstructor('Too many input arguments'); error(id,msg);
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