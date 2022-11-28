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
            % default constructor (empty object)
        
        elseif nargin==1
            % copy constructor
            if isa(varargin{1},'halfspace')
                obj = varargin{1};     
            else
                throw(CORAerror('CORA:wrongValue','first',"'halfspace' object or normal vector"));
            end
        %one input
        elseif nargin==2
            % set normal vector
            obj.c = varargin{1};
            % set distance to origin
            obj.d = varargin{2};
        elseif nargin > 2
            % too many input arguments
            throw(CORAerror('CORA:tooManyInputArgs',2));
        end
        
        % set parent object properties
        obj.dimension = length(obj.c);
    end
end

methods (Static = true)
    hs = generateRandom(varargin) % generates random halfspace
end


end

%------------- END OF CODE --------------