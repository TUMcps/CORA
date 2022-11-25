function E = and(obj,S,mode)
% and - overloads & operator to compute the intersection an ellipsoid and a
% set representation
%
% Syntax:  
%    [E] = and(obj,S)
%    [E] = and(obj,S,mode)
%
% Inputs:
%    obj            - Ellipsoid object
%    S              - set representation (or cell array thereof)
%    mode(optional) - approximation mode ('i':inner; 'o': outer)
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    E1= ellipsoid.generateRandom(2,false);
%    E2= ellipsoid.generateRandom(2,false);
%    E3 =ellipsoid.generateRandom(2,false);
%    E = E1 & E2;
%    E = and(E1,{E2,E3},'i');
%
% References: 
%   [1] A. Kurzhanski et al. "Ellipsoidal Toolbox Manual", 2006
%       https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%            
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  15-October-2019
%               15-March-2021
% Last revision:---

%------------- BEGIN CODE --------------
%% parsing and checking
if ~exist('mode','var')
    mode = 'o';
end
if ~any(mode==['i','o'])
   error('mode has to be either "i" (inner approx) or "o" (outer approx)');
end

% handle empty cells, S not a cell etc
S = prepareSetCellArray(S,obj);
if isempty(obj) || isempty(S) 
    E = ellipsoid;
    return;
end

N = length(S);
% if only center remains
if rank(obj)==0
    if ~ismethod(S{1},'in')
        error(['Ellipsoid contains only 1 point, but second argument type ',...
                'does not implement "in"!']);
    end
    if all(cellfun(@(s)in(s,obj.q),S))
        E = ellipsoid(zeros(dim(obj)),obj.q);
    else
        E = ellipsoid;
    end
    return;
end

%% different intersections

% ellipsoid and point
if isa(S{1},'double')
    S = cell2mat(reshape(S,[1,numel(S)]));
    % if not all points are equal, overall intersection is empty
    if ~all(all(withinTol(S,repmat(S(:,1),1,size(S,2)),obj.TOL))) ||...
        ~in(obj,S(:,1),'exact')
        E = ellipsoid;
    else
        E = ellipsoid(zeros(size(obj.Q)),S(:,1));
    end
    return;
end

% ellipsoid and conPolyZono
if isa(S{1},'conPolyZono')
    if strcmp(mode,'o')
        E = S{1} & obj;
        for i=2:N
            E = S{2} & E;
        end
    else
        error('inner approximation of Ellipsoid and conPolyZono not implemented');
    end
    return;
end

% ellipsoid and ellipsoid
if isa(S{1},'ellipsoid')
    if strcmp(mode,'o')
        E = andEllipsoidOA(obj,S{1});
        for i=2:N
            if isempty(E)
                break;
            end
            E = andEllipsoidOA(E,S{i});
        end
    else
        E = andEllipsoidIA(obj,S);
    end
    return;
end

% ellipsoid and conHyperplane
if isa(S{1},'conHyperplane')
    E = obj;
    for i=1:N
        if isHyperplane(S{i})
            E = andHyperplane(E,S{i});
        else
            E = and(E,mptPolytope(S{i}),mode);
        end
    end
    return;
    
end

% ellipsoid and mptPolytope
if isa(S{1},'mptPolytope')
    error('Not implemented yet');
end

% ellipsoid and halfspace
if isa(S{1},'halfspace')
    E = andHalfspace(obj,S{1},mode);
    for i=2:N
        E = andHalfspace(E,S,mode);
    end
    return;
end

error('Wrong type of input arguments');

%------------- END OF CODE --------------