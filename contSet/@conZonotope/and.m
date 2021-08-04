function res = and(obj,S)
% and - Computes the intersection of a constrained zonotope with
%             other set representations
%
% Syntax:  
%    res = and(obj,S)
%
% Inputs:
%    obj - constrained zonotope object
%    S - second set (supported objects: conZonotope, halfspace, 
%                    constrainedHyperlane)
%
% Outputs:
%    res - constrained zonotope object
%
% Example: 
%    % constrained zonotopes
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1];
%    b = 1;
%    cZono1 = conZonotope(Z,A,b);
%
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1];
%    b = 1;
%    cZono2 = conZonotope(Z,A,b);
%
%    % halfspace and constrained hyperplane
%    hs = halfspace([1,-2],1);
%    ch = conHyperplane([1,-2],1,[-2 -0.5;1 0],[-4.25;2.5]);
%
%    % compute intersection
%    res1 = cZono1 & cZono2;
%    res2 = cZono2 & hs;
%    res3 = cZono1 & ch;
%
%    % visualization
%    figure; hold on
%    plot(cZono1,[1,2],'r');
%    plot(cZono2,[1,2],'b');
%    plot(res1,[1,2],'g','Filled',true,'EdgeColor','none');
%    title('Constrained zonotope');
%
%    figure; hold on
%    xlim([-4,4]); ylim([-4,4]);
%    plot(hs,[1,2],'r','FaceAlpha',0.5);
%    plot(res2,[1,2],'g','Filled',true,'EdgeColor','none');
%    plot(cZono2,[1,2],'b');
%    title('halfspace');
%
%    figure; hold on
%    xlim([0,4]); ylim([-3,4]);
%    plot(ch,[1,2],'g');
%    plot(cZono1,[1,2],'r');
%    plot(res3,[1,2],'b','LineWidth',2);
%    title('Constrained hyperplane');    
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

% Author:        Dmitry Grebenyuk, Niklas Kochdumper
% Written:       13-May-2018
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: ---

%------------- BEGIN CODE --------------

    % get conZonotope object
    if ~isa(obj,'conZonotope')
       temp = obj;
       obj = S;
       S = temp;
    end
    
    % Add trivial constraint if the conZonotope object does not have
    % constraints (for easier implementation of the following operations)
    if isempty(obj.A)
       obj.A = zeros(1,size(obj.Z,2)-1);
       obj.b = 0;
    end
    
    % different cases depending on the set representation
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
        
        % delete all zero constraints and generators
        res = deleteZeros(res);
        

    elseif isa(S, 'halfspace')

        % Extract object properties C*x <= d of the halfspace
        C = S.c';
        d = S.d;

        G = obj.Z(:,2:end);
        c = obj.Z(:,1);

        % compute lower bound
        l = supportFunc(obj,C','lower');
        
        % Add additional constraints
        A = [obj.A, zeros(size(obj.A,1),1); C*G, 0.5*(l-d)];
        b = [obj.b; 0.5*(d+l)-C*c];
        G = [G,zeros(size(G,1),1)];

        res = conZonotope([c,G],A,b);


    elseif isa(S, 'conHyperplane')

        % calculate intersection between constrained zonotope and hyperplane
        C = (S.h.c)';
        d = S.h.d;

        G = obj.Z(:,2:end);
        c = obj.Z(:,1);

        A = [obj.A; C*G];
        b = [obj.b; d - C*c];

        res = conZonotope([c,G],A,b); 

        % loop over all constraints
        C = S.C;
        d = S.d;
        
        for i = 1:size(C,1)
            
           % construct halfspace
           hs = halfspace(C(i,:)',d(i));
           
           % check if set is fully contained in halfspace
           if ~in(hs,res)
              
               % intersect set with halfspace
               res = res & hs; 
           end
        end

    elseif isa(S,'mptPolytope')
        
        % get matrices for inequality constraints A*x <= b
        A = S.P.A;
        b = S.P.b;
        
        res = obj;
        
        % loop over all constraints
        for i = 1:size(A,1)
            
           % construct halfspace
           hs = halfspace(A(i,:)',b(i));
           
           % check if set is fully contained in halfspace
           if ~in(hs,res)
              
               % intersect set with halfspace
               res = res & hs; 
           end
        end        

    elseif isa(S,'zonotope') || isa(S,'interval') || isa(S,'zonoBundle')

        % convert to constrained zonotope
        res = obj & conZonotope(S);
        
    elseif isa(S,'levelSet') || isa(S,'conPolyZono')
        
        res = S & obj;
        
    else
        % throw error for given arguments
        error(noops(obj,S));
    end

end

%------------- END OF CODE --------------