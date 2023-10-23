function [obj,S] = findClassArg(obj1,obj2,classname)
% findClassArg - finds the obj with specified class name, or throws error
%    if neither object is of given class
%
% Syntax:
%    [obj,S] = findClassArg(obj1,obj2,classname)
%
% Inputs:
%    obj1,obj2 - set representation objects/numerical matrix or vector
%
% Outputs:
%    obj - class object
%    S   - other object
%
% Example:
%    obj1 = ellipsoid.generateRandom();
%    obj2 = zonotope.generateRandom();
%    [obj,S] = findClassArg(obj2,obj1,'ellipsoid');
%    % obj is now obj1, S is now obj2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Victor Gassmann
% Written:       04-July-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if isa(obj1,classname)
    obj = obj1;
    S = obj2;
elseif isa(obj2,classname)
    obj = obj2;
    S = obj1;
else
    % neither is of given class => throw error
    throw(CORAerror('CORA:wrongValue','first/second',classname));
end

% ------------------------------ END OF CODE ------------------------------
