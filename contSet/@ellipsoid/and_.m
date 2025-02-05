function E = and_(E,S,mode)
% and_ - overloads '&' operator to compute the intersection an ellipsoid
%    and another set representation
%
% Syntax:
%    E = and_(E,S)
%    E = and_(E,S,mode)
%
% Inputs:
%    E - ellipsoid object
%    S - contSet object, numeric, cell-array
%    mode(optional) - approximation mode ('inner','outer')
%
% Outputs:
%    E - ellipsoid object
%
% Example:
%    E1 = ellipsoid([3 -1; -1 1],[1;0]);
%    E2 = ellipsoid([5 1; 1 2],[1;-1]);
%    E3 = ellipsoid([0.6 -0.4; -0.4 2.2],[0.5;0]);
%    Eo = and(E1,{E2,E3},'outer');
%    Ei = and(E1,{E2,E3},'inner');
%    figure; hold on;
%    plot(E1); plot(E2); plot(E3);
%    plot(Eo); plot(Ei);
%
% References: 
%   [1] A. Kurzhanski et al. "Ellipsoidal Toolbox Manual", 2006
%       https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%            
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and

% Authors:       Victor Gassmann
% Written:       13-March-2019
% Last update:   15-October-2019
%                15-March-2021
%                04-July-2022 (VG, replaced cell arrays by class arrays)
% Last revision: 27-March-2023 (MW, rename and_)
%                28-September-2024 (MW, integrate precedence)

% ------------------------------ BEGIN CODE -------------------------------

% call function with lower precedence
if isa(S,'contSet') && S.precedence < E.precedence
    E = and_(S,E,mode);
    return
end

% read out dimension
n = dim(E);

% ellipsoid and point
if isnumeric(S) && iscolumn(S)
    if contains_(E,S,'exact',eps,0,false,false)
        E = ellipsoid(zeros(n),S);
    else
        E = ellipsoid.empty(n);
    end
    return;
end

% rewrite S as cell-array for unified handling
if ~iscell(S)
    S = {S};
end

% ellipsoid is only a point
if representsa_(E,'point',eps) 
    % if all sets contain that point, the intersection is that point
    if all(cellfun(@(S_i) contains_(S_i,E.q,'exact',eps,0,false,false), S, 'UniformOutput', true))
        E = ellipsoid(zeros(n),E.q);
    else
        E = ellipsoid.empty(n);
    end
    return;
end

% ellipsoid and polytope (including hyperplanes)
if all(cellfun(@(S_i) isa(S_i,'polytope'),S,'UniformOutput',true))
    for i=1:numel(S)
        if representsa_(S{i},'hyperplane',eps)
            % note: evaluation is exact
            E = priv_andHyperplane(E,S{i});
        else
            E = priv_andPolytope(E,S{i},mode);
        end
    end
    return;
end

% ellipsoid and ellipsoid
if all(cellfun(@(S_i) isa(S_i,'ellipsoid'),S,'UniformOutput',true))
    if strcmp(mode,'outer')
        for i=1:numel(S)
            E = priv_andEllipsoidOA(E,S{i});
            if representsa_(E,'emptySet',eps)
                return;
            end
        end
    else
        E = priv_andEllipsoidIA([{E}, S]);
    end
    return;
end

% throw error for remaining combinations
throw(CORAerror('CORA:noops',E,S));

% ------------------------------ END OF CODE ------------------------------
