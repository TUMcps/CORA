classdef halfspace
% halfspace class 
%
% Syntax:  
%    object constructor: Obj = halfspace(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    c - normal vector of the halfspace
%    d - distance to the origin (scalar)
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      06-June-2011
% Last update:  14-Aug-2019
%               02-May-2020 (add property validation)
% Last revision:---

%------------- BEGIN CODE --------------


properties (SetAccess = private, GetAccess = public)
    c (:,1) {mustBeNumeric,mustBeFinite} = [];
    d (1,1) {mustBeNumeric,mustBeFinite} = 0;
end
    
methods
    %class constructor
    function obj = halfspace(varargin)
        
        if nargin == 0
            % default constructor
        
        elseif nargin==1
            %copy constructor when halfspace function is called
            if isa(varargin{1},'halfspace')
                obj = varargin{1};     
            end
        %one input
        elseif nargin==2
            % set normal vector
            obj.c = varargin{1};
            % set distance to origin
            obj.d = varargin{2};
        end
    end
    
end

methods (Static = true)
    h = generateRandom(varargin) % generates random halfspace
end


end

%------------- END OF CODE --------------