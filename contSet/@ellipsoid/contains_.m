function res = contains_(E,S,type,varargin)
% contains_ - determines if an ellipsoid contains a set or a point
%
% Syntax:
%    res = contains_(E,S)
%    res = contains_(E,S,type)
%
% Inputs:
%    E - ellipsoid object 
%    S - contSet object or single point
%    type - mode of check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    E1 = ellipsoid([5 7;7 13],[1;2]);
%    E2 = ellipsoid(0.3*eye(2));
%    Z = zonotope([0 1 0;0 1 1]);
%
%    contains(E1,E2)
%    contains(E1,Z)
%
%    figure; hold on
%    plot(E1,[1,2],'b');
%    plot(E2,[1,2],'g');
%    plot(Z,[1,2],'r');
%
% References:
%    [1] Yildirim, E.A., 2006. On the minimum volume covering ellipsoid of
%        of ellipsoids. SIAM Journal on Optimization, 17(3), pp.621-641.     
%    [2] SDPT3: url: http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains

% Authors:       Victor Gassmann, Niklas Kochdumper, Adrian Kulmburg
% Written:       15-October-2019 
% Last update:   21-November-2019 (NK, extend to other sets)
%                09-March-2021 (included tolerance for q comparison)
%                17-March-2021 (error handling)
%                06-July-2021 (AK, merged containsPoint to in)
%                04-July-2022 (VG, handle class array cases)
%                25-November-2022 (MW, rename 'contains')
%                25-April-2023 (VG, add method for capsule)
% Last revision: 27-March-2023 (MW, rename contains_)

% ------------------------------ BEGIN CODE -------------------------------

% containment check
if isnumeric(S)
    res = aux_containsPoint(E,S);

elseif isa(S,'ellipsoid')
    res = containsEllipsoid(E,S);

elseif isa(S,'capsule')
    % check if balls at both ends of capsule are contained
    E1 = ellipsoid(S.r^2*eye(dim(S)), S.c+S.g);
    E2 = ellipsoid(S.r^2*eye(dim(S)), S.c-S.g);
    res = containsEllipsoid(E,E1) && containsEllipsoid(E,E2);

elseif isa(S,'zonotope') && strcmp(type,'approx')
    S = ellipsoid(S,'outer:norm_bnd');
    res = containsEllipsoid(E,S);

else
    if strcmp(type,'exact')
        if ismethod(S,'vertices_')
            % check if all vertices of the set are contained
            res = all(contains_(E,vertices(S),type,0));
        else
            throw(CORAerror('CORA:noExactAlg',E,S));
        end   
    else
        if ismethod(S,'zonotope')
            % zonotope over-approximation
            S = zonotope(S); 
            res = contains_(E,S,type,0);
        else
            throw(CORAerror('CORA:noops',E,S));
        end
    end

end

end


% Auxiliary functions -----------------------------------------------------

function [res,val] = aux_containsPoint(E,Y)
% gives an array of boolean values indicating whether points Y are
% contained in the ellipsoid; additionally, the relative distance to the
% center of the ellipsoid is returned (<= 1 if contained, =1 on boundary,
% otherwise Inf)

% number of points
N = size(Y,2);

% init results
res = false(1,N);
val = inf(1,N);
ind_rem_eq = true(1,N);

% degenerate case
if ~isFullDim(E)

    % subspace dimension
    nt = rank(E);

    % all-zero Q matrix -> ellipsoid is just a point
    if nt == 0
        res = all(withinTol(repmat(E.q,1,N),Y,E.TOL),1);
        val(res) = 0;
        return;
    end

    % project via svd
    [T,~,~] = svd(E.Q);
    E = T'*E;
    Y = T'*Y;

    % save remainder
    x_rem = E.q(nt+1:end);
    Y_rem = Y(nt+1:end,:);

    % check whether x_rem==Y_rem (those that do not fulfill are already
    % not contained)

    % indices of res which might be contained
    ind_rem_eq = all(withinTol(repmat(x_rem,1,size(Y_rem,2)),Y_rem,E.TOL),1);

    % if only center remains
    if rank(E)==0
        res(ind_rem_eq) = true;
        val(ind_rem_eq) = 1;
        return;
    end

    % project so that E is no longer degenerate
    rankE = rank(E);
    E = project(E,1:rankE);
    Y = Y(1:rankE,:);
end

% convert mask to indices
tmp = 1:N;
ii_eq_rem = tmp(ind_rem_eq);

% now, E is full-dimensional
for i=ii_eq_rem
    % simply check using ellipsoid equation
    val_i = (Y(:,i)-E.q)'*(E.Q\(Y(:,i)-E.q));
    tmp_E = 1 + E.TOL;
    res(i) = val_i < tmp_E | withinTol(val_i,tmp_E);
    val(i) = res(i)*val_i;
end

end

% ------------------------------ END OF CODE ------------------------------
