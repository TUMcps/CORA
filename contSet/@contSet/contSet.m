classdef contSet
% contSet - Object and Copy Constructor 
%
% Syntax:  
%    object constructor: Obj = class_name(varargin)
%    copy constructor: Obj = otherObj
%
% Inputs:
%    input1 - dimension: int
%
% Outputs:
%    Obj - Generated Object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:        Matthias Althoff, Mark Wetzlinger
% Written:       02-May-2007 
% Last update:   04-May-2020 (MW, transition to classdef)
%                01-June-2022 (MW, add CORAerror)
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    dimension (1,1) {mustBeInteger} = 0;
end

methods

    function obj = contSet(varargin)
        
        % (default constructor)
        if nargin == 0
            % values in properties already initialized
            
        % If 1 argument is passed
        elseif nargin == 1
            % (copy constructor)
            if isa(varargin{1}, 'contSet')
                obj = varargin{1};
            else
                %List elements of the class
                obj.dimension = varargin{1};
            end

        % Else if not enough or too many inputs are passed    
        else
            throw(CORAerror('CORA:tooManyInputArgs',1));
        end
    end
end
end

%------------- END OF CODE --------------