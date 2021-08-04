function res = in(obj1,obj2,mode)
% in - checks whether obj2 is contained in obj1
%
% Syntax:  
%    res = in(obj1,obj2)
%    res = in(obj,obj2,mode)
%
% Inputs:
%    obj1 - ellipsoid object 
%    obj2 - contSet object 
%    mode - mode of check ('exact' or 'approx')
%
% Outputs:
%    res - result of containment check
%
% Example: 
%    E1 = ellipsoid([5 7;7 13],[1;2]);
%    E2 = ellipsoid(0.3*eye(2));
%    zono = zonotope([0 1 0;0 1 1]);
%
%    in(E1,E2)
%    in(E1,zono)
%
%    figure
%    hold on
%    plot(E1,[1,2],'b');
%    plot(E2,[1,2],'g');
%
%    figure
%    hold on
%    plot(E1,[1,2],'b');
%    plot(zono,[1,2],'r');
%
% References:
%            [1] Yildirim, E.A., 2006. On the minimum volume covering 
%                ellipsoid of ellipsoids. SIAM Journal on Optimization, 
%                17(3), pp.621-641.     
%            [2] SDPT3: url: http://www.math.nus.edu.sg/~mattohkc/sdpt3.html
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann, Niklas Kochdumper, Adrian Kulmburg
% Written:      15-October-2019 
% Last update:  21-November-2019 (NK, extend to other sets)
%               09-March-2021 (included tolerance for q comparison)
%               17-March-2021 (error handling)
%               06-July-2021 (AK, merged containsPoint to in)
% Last revision:---

%------------- BEGIN CODE --------------
if ~exist('mode','var')
    mode = 'exact';
end
if isa(obj2,'double')

    for i = 1:size(obj2,2)
        res = containsPoint(obj1,obj2(:,i)); 
        if ~res
            return;
        end
    end
elseif isa(obj2,'ellipsoid')
    res = inEllipsoid(obj1,obj2);
elseif isa(obj2,'zonotope') && strcmp(mode,'approx')
    obj2 = ellipsoid(obj2,'o:norm:bnd');
    res = inEllipsoid(obj1,obj2);
else
    if strcmp(mode,'exact')
        if ismethod(obj2,'vertices')
            % check if all vertices of the set are contained
            res = in(obj1,vertices(obj2));
        else
            error('Mode "exact" not implemented for second argument type!');
        end   
    else
        if ismethod(obj2,'zonotope')
            % zonotope over-approximation
            obj2 = zonotope(obj2); 
            res = in(obj1,obj2);
        else
            error('Not implemented for second argument type!');
        end
    end
end

function [B,Val] = containsPoint(E,Y)
% containsPoint - gives an array of boolean values indiciating whether
%    points Y are contained in the ellipsoid
%
% Syntax:  
%    [B,Val] = containsPoint(E,Y) gives an array of boolean values indiciating
%     whether points Y are contained in the ellipsoid
%
% Inputs:
%    E - ellipsoids object
%    Y - Points
%
% Outputs:
%    B - boolean values indiciating whether
%        points Y are contained in the ellipsoid
%    Val-if contained, Value indicates the relative distance to the center of E: Val<=1
%    <=> contained (=1: on boundary); otherwise: inf
%
% Example: 
%    t = linspace(0,2*pi,1000);
%    Y = [cos(t);sin(t)];
%    E = ellipsoid([1,0;0,1/2],[1;1]);
%    B = containsPoint(E,Y);
%

n = dim(E);
if ~isa(Y,'double')
    error('Second argument must be a double matrix');
end
[m,N] = size(Y);
if m~=n
    error('First dimension of second input does not match ellipsoid dimension!');
end


B = false(1,N);
Val = inf(1,N);
ind_rem_eq = true(1,N);
if E.isdegenerate
    nt = rank(E);
    % if all-zero Q matrix
    if nt==0
        B = all(withinTol(repmat(E.q,1,N),Y,E.TOL),1);
        Val(B) = 0;
        return;
    end
    [T,~,~] = svd(E.Q);
    E = T'*E;
    Y = T'*Y;
    % save remainder
    x_rem = E.q(nt+1:end);
    Y_rem = Y(nt+1:end,:);
    % check whether x_rem==Y_rem (those that do not fullfill that are
    % already not contained)
    % indices of B which might be contained
    ind_rem_eq = all(withinTol(repmat(x_rem,1,size(Y_rem,2)),Y_rem,E.TOL));
    % if only center remains
    if rank(E)==0
        B(ind_rem_eq) = true;
        Val(ind_rem_eq) = 1;
        return;
    end
    % project so that E is no longer degenerate
    E = project(E,1:E.rank);
    Y = Y(1:E.rank,:);
end
% convert mask to indices
tmp = 1:N;
ii_eq_rem = tmp(ind_rem_eq);
% now, E is fulldimensional
for i=ii_eq_rem
    % simply check using ellipsoid equation
    val_i = (Y(:,i)-E.q)'*inv(E.Q)*(Y(:,i)-E.q);
    B(i) = val_i <= 1+E.TOL;
    Val(i) = B(i)*val_i;
end
%------------- END OF CODE --------------