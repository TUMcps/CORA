function res = in(obj,S)
% in - set containment operator for Signal Temporal Logic
%
% Syntax:  
%    res = in(obj,S)
%
% Inputs:
%    obj - logic formula (class stl)
%    S - set (class interval, zonotope, mptPolytope, conZonotope,
%               zonoBundle, halfspace, or polygon)
%
% Outputs:
%    res - resulting stl formula (class stl)
%
% Example: 
%    x = stl('x',2)
%    pgon = polygon.generateRandom();
%    eq = finally(in(x,pgon),interval(0.1,0.2))
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: stl

% Author:       Niklas Kochdumper
% Written:      09-November-2022 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% check input arguments
if ~isa(obj,'stl') || ~strcmp(obj.type,'variable')
    throw(CORAerror('CORA:notSupported',...
                  'This operation is not supported for stl objects!'));
end

if dim(S) ~= length(obj.variables)
    throw(CORAerror('CORA:wrongInput','second', ...
                  'dimensions of set and stl object have to match"!'));
end

% conversions for different types of sets
if isa(S,'polygon')

    list = splitIntoConvexSets(S);

    res = [];
    for i = 1:length(list)
        res = res | polytope2stl(obj,mptPolytope(list{i}));
    end

elseif isa(S,'interval') || isa(S,'zonotope') || ...
       isa(S,'mptPolytope') || isa(S,'conZonotope') || ...
       isa(S,'zonoBundle') || isa(S,'halfspace')

    res = polytope2stl(obj,mptPolytope(S));

else
    throw(CORAerror('CORA:notSupported',...
         'This operation is not supported for this type of set representation!'));
end

%------------- END OF CODE --------------