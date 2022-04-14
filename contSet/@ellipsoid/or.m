function E = or(obj,S,mode)
% or - overloads | operator to compute inner/outer approximation of the
% union 
%
% Syntax:  
%    E = or(obj,S)
%    E = or(obj,S,mode)
%
% Inputs:
%    obj                - ellipsoid object
%      S                - set (or cell array)
%    mode (optional)    - 'i' (inner)/ 'o' (outer)
%
% Outputs:
%    E - resulting ellipsoid
%
% Example: 
%    E1 = ellipsoid.generateRandom(2,false);
%    E2 = ellipsoid.generateRandom(2,false); 
%    E = E1 | E2;
%
% References:
%    Convex Optimization; Boyd
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      09-March-2021 
% Last update:  15-March-2021
% Last revision:---

%------------- BEGIN CODE --------------
%% parsing & checking
if ~exist('mode','var')
    mode = 'o';
end
if ~any(mode==['i','o'])
   error('mode has to be either "i" (inner approx) or "o" (outer approx)');
end

if ~isa(obj,'ellipsoid')
    error('First argument is not an ellipsoid!');
end

S = prepareSetCellArray(S,obj);
if isempty(S)
    E = ellipsoid;
    return;
end

%% different unions
if isa(S{1},'double')
    V = cell2mat(reshape(S,[1,numel(S)]));
    if strcmp(mode,'o')
        E_p = ellipsoid.enclosePoints(V,'min-vol');
        E = orEllipsoidOA(obj,{E_p});
        return;
    else
        P = mptPolytope.enclosePoints(V);
        E = ellipsoid(P,'i');
        return;
    end
end

if isa(S{1},'ellipsoid')
    if strcmp(mode','o')
        E = orEllipsoidOA(obj,S);
        return;
    else
        S_ = [{obj};S];
        [~,ii] = max(cellfun(@(ss)volume(ss),S_));
        E = S_{ii};
        return;
    end
end

if isa(S{1},'mptPolytope')
    error('Not implemented yet!');
end

error('Type of S not supported');
%------------- END OF CODE --------------