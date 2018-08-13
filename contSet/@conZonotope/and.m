function res = and(obj,S)
% intersect - Computes the intersection of a constrained zonotope with
%             other set representations
%
% Syntax:  
%    res = intersect(obj,Z)
%
% Inputs:
%    obj - constrained zonotope object
%    Z - second set (supported objects: conZonotope, halfspace, 
%                    constrainedHyperlane)
%
% Outputs:
%    res - constrained zonotope object
%
% Example: 
%    % constrained zonotope 1
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1];
%    b = 1;
%    cZono1 = conZonotope(Z,A,b);
%
%    % constrained zonotope 2
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1];
%    b = 1;
%    cZono2 = conZonotope(Z,A,b);
%
%    % hyperplane
%    C = [1 -2];
%    d = 1;
%    hp = halfspace(C,d);
%
%    % constrained hyperplane
%    Ch = [-2 -0.5;1 0];
%    dh = [-4.25;2.5];
%    ch = constrainedHyperplane(hp,Ch,dh);
%
%    % intersection between two constrained zonotopes
%    intZono = cZono1 & cZono2;
%    figure
%    plot(cZono1,[1,2],'r');
%    hold on
%    plot(cZono2,[1,2],'b');
%    plotFilled(intZono,[1,2],'g');
%    title('Constrained zonotope');
%
%    % insection with hyperplane
%    intZono = cZono2 & hp;
%    x = -4:0.1:4;
%    y = (d-C(1)*x)./C(2);
%    figure
%    hold on
%    plot(x,y,'g');
%    plot(cZono2,[1,2],'r');
%    plot(intZono,[1,2],'b');
%    title('hyperplane');
%
%    % intersection with constrained hyperplane
%    intZono = cZono1 & ch;
%    figure
%    hold on
%    poly = mptPolytope([Ch;0 1],[dh;4]);
%    plotFilled(poly,[1,2],'m','EdgeColor','none','FaceAlpha',0.5);
%    plot(x,y,'g');
%    plot(cZono1,[1,2],'r');
%    plot(intZono,[1,2],'b','LineWidth',2);
%    title('Constrained hyperplane');
%       
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Author: Dmitry Grebenyuk, Niklas Kochdumper
% Written: 13-May-2018
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

if ~isempty(S)
    
    % Add trivial constraint if the conZonotope object does not have
    % constraints (for easier implementation of the following operations)
    if isempty(obj.A)
       obj.A = zeros(1,size(obj.Z,2)-1);
       obj.b = 0;
    end
    
    
    if isa(S, 'conZonotope') 
        
        if isempty(S.A)
           S.A = zeros(1,size(S.Z,2)-1);
           S.b = 0;
        end
        
        % Calculate intersection according to equation (13) at Proposition 1 in
        % reference paper [1]
        Z = [obj.Z, zeros(size(S.Z)-[0,1])];
        A = blkdiag(obj.A,S.A);
        A = [A; obj.Z(:,2:end), -S.Z(:,2:end)];
        b = [obj.b; S.b; S.Z(:,1) - obj.Z(:,1)];

        res = conZonotope(Z,A,b);


    elseif isa(S, 'zonotope')
        
        % Calculate intersection according to equation (13) at Proposition 1 in
        % reference paper [1]
        Z = [obj.Z, zeros(size(S.Z)-[0,1])];
        A = [obj.A zeros(size(obj.A,1),size(S.Z,2)-1); obj.Z(:,2:end), -S.Z(:,2:end)];
        b = [obj.b; S.Z(:,1) - obj.Z(:,1)];

        res = conZonotope(Z,A,b);


    elseif isa(S, 'halfspace')

        % Extract object properties C*x = d of the hyperplane
        C = get(S,'c')';
        d = get(S,'d');

        G = obj.Z(:,2:end);
        c = obj.Z(:,1);

        % Add additional constraints
        A = [obj.A; C*G];
        b = [obj.b; d - C*c];

        res = conZonotope([c,G],A,b);


    elseif isa(S, 'constrainedHyperplane')

        % Calculate intersection between constrained zonotope and hyperplane
        res = obj & S.h; 

        % Check if the original unconstrained zonotope violates the constraints
        % (fast test to see if the computational expensive instersection with
        % the constraits has to be performed)
        zono = zonotope(obj.Z);

        if ~constrSat(zono, S.C, S.d)    % constraints have to be considered

           % calculate a bounding box for the unconstrained zonotope
           int = interval(zono);

           % calculate a polytope that represents the intersection of the
           % bounding box with the feasible region of the constraints
           % (polytope has to be closed for numerical stability)
           n = length(int);
           A = [-eye(n);eye(n);S.C];
           b = [-infimum(int);supremum(int);S.d];

           poly = mptPolytope(A,b);

           % intersect the polytope with the part on the hyperplane
           cZono = conZonotope(poly);
           res = res & cZono;     
        end


    elseif isa(S,'mptPolytope')

        % convert polytope to constrained zonotope
        res = obj & conZonotope(S);

    elseif isa(S,'interval')
        
        % convert interval to  zonotope
        res = obj & zonotope(S);
        
    else
        error('Wrong input')
    end
else
    res = conZonotope([],[],[]);
end
end

%------------- END OF CODE --------------