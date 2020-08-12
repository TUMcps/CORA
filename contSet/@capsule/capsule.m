classdef capsule
% capsule - Object and Copy Constructor 
%
% Syntax:  
%    obj = capsule(c)
%    obj = capsule(r)
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
%    obj - Generated Object
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
% See also: interval,  polytope

% Author:       Matthias Althoff
% Written:      04-March-2019
% Last update:  02-May-2020 (MW, add property validation)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    c (:,1) {mustBeNumeric,mustBeFinite} = []; % center
    g = []; % generator
    r (1,1) {mustBeNonnegative,mustBeFinite} = 0; % radius
    contSet = [];
end
   
methods

    function Obj = capsule(varargin)
        
        % default constructor
        if nargin == 0
            
            Obj.contSet = contSet();
            
        else
            
            % If 1 argument is passed
            if nargin == 1
                if isa(varargin{1},'capsule')
                    % copy constructor
                    Obj = varargin{1};
                else
                    if isvector(varargin{1})
                        % set center
                        Obj.c = varargin{1}; 
                    elseif isscalar(varargin{1})
                        % set radius
                        Obj.r = varargin{1}; 
                    end
                end
            elseif nargin == 2
                if isvector(varargin{1}) && isvector(varargin{2})
                    % set center and generator
                    Obj.c = varargin{1}; 
                    Obj.g = varargin{2};
                elseif isvector(varargin{1}) && isscalar(varargin{2})
                    % set center and radius
                    Obj.c = varargin{1}; 
                    Obj.r = varargin{2};
                end
            elseif nargin == 3
                % set all values
                Obj.c = varargin{1}; 
                Obj.g = varargin{2};
                Obj.r = varargin{3};
            end

            %Generate parent object
            if ~isempty(varargin{1}) && ~isa(varargin{1}, 'capsule') && isvector(varargin{1})
                Obj.contSet = contSet(length(varargin{1}(:,1)));
            else
                Obj.contSet = contSet();
            end
        
        end
        
    end
end

methods (Static = true)
    C = generateRandom(varargin) % generates random capsule
end


end
%------------- END OF CODE --------------