function res = and_(cZ,S,varargin)
% and_ - Computes the intersection of a constrained zonotope with
%    other set representations
%
% Syntax:
%    res = and_(cZ,S)
%
% Inputs:
%    cZ - conZonotope object
%    S - contSet object
%
% Outputs:
%    res - conZonotope object
%
% Example: 
%    % constrained zonotopes
%    Z = [0 3 0 1;0 0 2 1];
%    A = [1 0 1]; b = 1;
%    cZ1 = conZonotope(Z,A,b);
%    Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
%    A = [1 1 1]; b = 1;
%    cZ2 = conZonotope(Z,A,b);
%
%    % halfspace and constrained hyperplane
%    hs = halfspace([1,-2],1);
%    hyp = conHyperplane([1,-2],1,[-2 -0.5;1 0],[-4.25;2.5]);
%
%    % compute intersection
%    res1 = cZ1 & cZ2;
%    res2 = cZ2 & hs;
%    res3 = cZ1 & hyp;
%
%    % visualization
%    figure; hold on;
%    plot(cZ1,[1,2],'r');
%    plot(cZ2,[1,2],'b');
%    plot(res1,[1,2],'FaceColor','g');
%    title('Constrained zonotope');
%
%    figure; hold on; xlim([-4,4]); ylim([-4,4]);
%    plot(hs,[1,2],'r','FaceAlpha',0.5);
%    plot(res2,[1,2],'FaceColor','g');
%    plot(cZ2,[1,2],'b');
%    title('halfspace');
%
%    figure; hold on; xlim([0,4]); ylim([-3,4]);
%    plot(hyp,[1,2],'g');
%    plot(cZ1,[1,2],'r');
%    plot(res3,[1,2],'b','LineWidth',2);
%    title('Constrained hyperplane'); 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and
%
% References: 
%   [1] J. Scott et al. "Constrained zonotope: A new tool for set-based
%       estimation and fault detection"

% Authors:       Dmitry Grebenyuk, Niklas Kochdumper
% Written:       13-May-2018
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: 27-March-2023 (MW, rename and_)

% ------------------------------ BEGIN CODE -------------------------------

% Add trivial constraint if the conZonotope object does not have
% constraints (for easier implementation of the following operations)
if isempty(cZ.A)
    cZ.A = zeros(1,size(cZ.G,2));
    cZ.b = 0;
end

% different cases depending on the set representation
if isa(S, 'conZonotope') 
    
    if isempty(S.A)
       S.A = zeros(1,size(S.G,2));
       S.b = 0;
    end
    
    % Calculate intersection according to equation (13) at Proposition 1 in
    % reference paper [1]
    Z = [cZ.c, cZ.G, zeros(size(S.G))];
    A = blkdiag(cZ.A,S.A);
    A = [A; cZ.G, -S.G];
    b = [cZ.b; S.b; S.c - cZ.c];

    res = conZonotope(Z,A,b);
    
    % delete all zero constraints and generators
    res = compact_(res,'zeros',eps);
    

elseif isa(S, 'halfspace')

    % Extract object properties C*x <= d of the halfspace
    C = S.c';
    d = S.d;

    G = cZ.G;
    c = cZ.c;

    % compute lower bound
    l = supportFunc_(cZ,C','lower');
    
    % Add additional constraints
    A = [cZ.A, zeros(size(cZ.A,1),1); C*G, 0.5*(l-d)];
    b = [cZ.b; 0.5*(d+l)-C*c];
    G = [G,zeros(size(G,1),1)];

    res = conZonotope([c,G],A,b);


elseif isa(S, 'conHyperplane')

    % calculate intersection between constrained zonotope and hyperplane
    C = S.a;
    d = S.b;

    G = cZ.G;
    c = cZ.c;

    A = [cZ.A; C*G];
    b = [cZ.b; d - C*c];

    res = conZonotope([c,G],A,b); 

    % loop over all constraints
    C = S.C;
    d = S.d;
    
    for i = 1:size(C,1)
        
       % construct halfspace
       hs = halfspace(C(i,:)',d(i));
       
       % check if set is fully contained in halfspace
       if ~contains_(hs,res)
          
           % intersect set with halfspace
           res = and_(res,hs,'exact'); 
       end
    end

elseif isa(S,'polytope')
    
    % get matrices for inequality constraints A*x <= b
    A = S.A;
    b = S.b;
    
    res = cZ;
    
    % loop over all constraints
    for i = 1:size(A,1)
        
       % construct halfspace
       hs = halfspace(A(i,:)',b(i));
       
       % check if set is fully contained in halfspace
       if ~contains_(hs,res)
          
           % intersect set with halfspace
           res = and_(res,hs,'exact');
       end
    end        

elseif isa(S,'zonotope') || isa(S,'interval') || isa(S,'zonoBundle')

    % convert to constrained zonotope
    res = and_(cZ,conZonotope(S),'exact');
    
elseif isa(S,'levelSet') || isa(S,'conPolyZono')
    
    res = and_(S,cZ,'exact');
    
else
    % throw error for given arguments
    throw(CORAerror('CORA:noops',cZ,S));
end

% ------------------------------ END OF CODE ------------------------------
