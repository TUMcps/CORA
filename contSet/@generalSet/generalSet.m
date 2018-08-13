classdef generalSet
% generalSet class 
%
% Syntax:  
%    object constructor: obj = generalSet(varargin)
%    copy constructor: obj = otherObj
%
% Inputs:
%    input1 - function handle
%    input2 - variable set
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval,  polytope

% Author:       Matthias Althoff
% Written:      04-November-2010
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    fct_handle = [];
    variable_set = [];
    segmentation = [];
    samples = [];
end
    
methods
    %class constructor
    function obj = generalSet(handle,contSet,seg)
        obj.fct_handle = handle;
        obj.variable_set = contSet;
        obj.segmentation = seg;
    end
         
    %methods in seperate files 
    Z = parallelotope(obj)
    
    %display functions
    plot(varargin)
    display(obj)
end
end

%------------- END OF CODE --------------