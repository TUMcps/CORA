function E = and(E,S,varargin)
% and - overloads '&' operator to compute the intersection an ellipsoid and
%    another set representation
%
% Syntax:  
%    E = and(E,S)
%    E = and(E,S,mode)
%
% Inputs:
%    E              - ellipsoid object
%    S              - set representation (array)
%    mode(optional) - approximation mode ('inner','outer')
%
% Outputs:
%    E - ellipsoid object
%
% Example: 
%    E1 = ellipsoid.generateRandom('Dimension',2);
%    E2 = ellipsoid.generateRandom('Dimension',2);
%    E3 = ellipsoid.generateRandom('Dimension',2);
%    Eo = and(E1,[E2,E3],'outer');
%    Ei = and(E1,[E2,E3],'inner');
%    figure; hold on;
%    plot(E1); plot(E2);plot(E3);
%    plot(Eo,[1,2],'r');
%    plot(Eo,[1,2],'b');
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
%               04-July-2022 (VG: replaced cell arrays by class arrays)
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_and('ellipsoid',E,S,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    E = vars{1}; return
else
    % potential re-ordering
    E = vars{1}; S = vars{2}; mode = vars{3};
end


N = length(S);
% if only center remains
if rank(E)==0 && isa(S,'contSet')
    % if double, already taken care if in 93
    if ~ismethod(S,'in')
        throw(CORAerror('CORA:noops',E,S{1}));
    end
    if all(arrayfun(@(ii)contains(S(ii),E.q),1:N))
        E = ellipsoid(zeros(dim(E)),E.q);
    else
        E = ellipsoid;
    end
    return;
end

%% different intersections

% ellipsoid and point
if isa(S,'double')
    % if not all points are equal, overall intersection is empty
    if ~all(all(withinTol(S,repmat(S(:,1),1,size(S,2)),E.TOL))) || ...
        ~contains(E,S(:,1),'exact')
        E = ellipsoid;
    else
        E = ellipsoid(zeros(size(E.Q)),S(:,1));
    end
    return;
end

% ellipsoid and conPolyZono
if isa(S,'conPolyZono')
    if strcmp(mode,'outer')
        E = S(1) & E;
        for i=2:N
            E = S(i) & E;
        end
    else
        throw(CORAerror('CORA:noops',E,S));
    end
    return;
end

% ellipsoid and ellipsoid
if isa(S,'ellipsoid')
    if strcmp(mode,'outer')
        E = andEllipsoidOA(E,S(1));
        for i=2:N
            if isempty(E)
                break;
            end
            E = andEllipsoidOA(E,S(i));
        end
    else
        E = andEllipsoidIA(E,S);
    end
    return;
end

% ellipsoid and conHyperplane
if isa(S,'conHyperplane')
    for i=1:N
        if isHyperplane(S(i))
            E = andHyperplane(E,S(i));
        else
            E = and(E,mptPolytope(S(i)),mode);
        end
    end
    return;
end

% ellipsoid and mptPolytope
if isa(S,'mptPolytope')
    E = andMptPolytope(E,S(1),mode);
    for i=2:N
        E = andMptPolytope(E,S(i),mode);
    end
    return;
end

% ellipsoid and halfspace
if isa(S,'halfspace')
    E = andHalfspace(E,S(1),mode);
    for i=2:N
        E = andHalfspace(E,S(i),mode);
    end
    return;
end

% throw error for remaining combinations
throw(CORAerror('CORA:noops',E,S));

%------------- END OF CODE --------------