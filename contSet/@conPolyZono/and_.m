function res = and_(cPZ,S,varargin)
% and_ - Computes the intersection of a constrained polynomial zonotope and
%    other set representations
%
% Syntax:
%    res = and_(cPZ,S)
%
% Inputs:
%    cPZ - conPolyZono object
%    S - contSet object
%
% Outputs:
%    res - conPolyZono object
%
% Example: 
%    % intersection with halfspace
%    cPZ = conPolyZono([0;0],[2 0 2;0 2 2],[1 0 3;0 1 1]);
%    hs = halfspace([1 1],0);
%
%    res1 = cPZ & hs;
%
%    figure; hold on;
%    xlim([-5,5]); ylim([-5,5]);
%    plot(hs,[1,2],'r');
%    plot(res1,[1,2],'FaceColor','b','Splits',10);
%    plot(cPZ,[1,2],'g','LineWidth',2);
%
%    % intersection with constrained hyperplane
%    cPZ = conPolyZono([0;0],[2 0 2;0 2 2],[1 0 3;0 1 1]);
%    hyp = conHyperplane([1 1],0,[-1 0],1);
%
%    res2 = cPZ & hyp;
%
%    figure; hold on; xlim([-5,5]); ylim([-5,5]);
%    plot(hyp,[1,2],'r');
%    plot(res2,[1,2],'FaceColor','b','Splits',10);
%    plot(cPZ,[1,2],'g','LineWidth',2);
%
%    % intersection of two constrained polynomial zonotopes
%    cPZ1 = conPolyZono([0;0],[2 0 2;0 2 2],[1 0 3;0 1 1]);
%    cPZ2 = conPolyZono([-1;1],[2 0 2;0 4 -4],[1 0 2;0 1 1]);
%
%    res3 = cPZ1 & cPZ2;
%
%    figure; hold on;
%    plot(res3,[1,2],'FaceColor',[0 .5 0],'Splits',10);
%    plot(cPZ1,[1,2],'r','LineWidth',1.5);
%    plot(cPZ2,[1,2],'b','LineWidth',1.5);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and, conZonotope/and_, interval/and_

% Authors:       Niklas Kochdumper
% Written:       07-November-2018
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename and_)

% ------------------------------ BEGIN CODE -------------------------------

% different cases for different set representations
if isa(S,'conPolyZono') || isa(S,'polytope') || ...
   isa(S,'interval') || isa(S,'zonotope') || ...
   isa(S,'zonoBundle') || isa(S,'conZonotope') || ...
   isa(S,'ellipsoid') || isa(S,'capsule') || ...
   isa(S,'polyZonotope') || isa(S,'taylm')

    % convert to constrained polynomial zonotope
    S = conPolyZono(S);
    
    % compute constraint part of the resulting set
    cPZ.b = [cPZ.b; S.b; S.c - cPZ.c];
    A = [cPZ.G  -S.G];
    EC = blkdiag(cPZ.E,S.E);
    
    if isempty(cPZ.A) && isempty(S.A)         
        cPZ.A = A;
        cPZ.EC = EC;
    elseif isempty(cPZ.A)
        cPZ.A = blkdiag(S.A,A);
        O = zeros(size(cPZ.E,1),size(S.EC,2));
        cPZ.EC = [[O;S.EC],EC];
    elseif isempty(S.A)
        cPZ.A = blkdiag(cPZ.A,A);
        O = zeros(size(S.E,1),size(cPZ.EC,2));
        cPZ.EC = [[cPZ.EC;O],EC];
    else           
        cPZ.A = blkdiag(cPZ.A,S.A,A);
        O1 = zeros(size(S.E,1),size(cPZ.EC,2));
        O2 = zeros(size(cPZ.E,1),size(S.EC,2)); 
        cPZ.EC = [[cPZ.EC;O1],[O2;S.EC],EC];
    end
    
    % compute state part of the resulting set
    O = zeros(size(S.E,1),size(cPZ.E,2));
    cPZ.E = [cPZ.E;O];
    
    % compute independent part of the resulting set
    if ~isempty(cPZ.GI)
        if ~isempty(S.GI)
           n = dim(cPZ); cen = zeros(n,1);
           Z1 = zonotope(cen,cPZ.GI);
           Z2 = zonotope(cen,S.GI);
           Z = enclose(Z1,Z2);
           cPZ.GI = Z.G;
        end
    else
        cPZ.GI = S.GI;
    end
    
    % update identifier vector
    cPZ.id = [cPZ.id; (max(cPZ.id)+1:max(cPZ.id)+length(S.id))'];
    
    % remove redundant monomials
    res = compact_(cPZ,'all',eps);

    
elseif isa(S,'halfspace')

    % compute lower bound
    l = supportFunc_(cPZ,S.c,'lower','interval',[]);
    
    % add additional constraints
    A = [S.c' * cPZ.G, -0.5*(S.d-l)];
    b = 0.5*(S.d+l) - S.c'*cPZ.c;
    EC = blkdiag(cPZ.E,1);
    
    cPZ.E = [cPZ.E; zeros(1,size(cPZ.E,2))];
    cPZ.id = [cPZ.id;max(cPZ.id) + 1];
    
    if isempty(cPZ.A)
       cPZ.A = A; cPZ.b = b; cPZ.EC = EC;
    else
       cPZ.A = blkdiag(cPZ.A, A);
       cPZ.b = [cPZ.b; b];
       cPZ.EC = [[cPZ.EC;zeros(1,size(cPZ.EC,2))], EC];
    end
    
    % remove redundant monomials
    res = compact_(cPZ,'all',eps);

    
elseif isa(S, 'conHyperplane')

    % calculate intersection between conPolyZono object and hyperplane
    % defined as c * x = d
    cPZ.b = [cPZ.b; S.b - S.a * cPZ.c];
    cPZ.EC = [cPZ.EC cPZ.E];
    
    if isempty(cPZ.A)
       cPZ.A = S.a * cPZ.G; 
    else
       cPZ.A = blkdiag(cPZ.A, S.a * cPZ.G); 
    end
    
    res = cPZ;

    % check if the zonotope over-approximation violates the constraints
    % (fast test to see if the computational expensive instersection 
    % with the constraits has to be performed)
    Z = zonotope(res);

    for i = 1:size(S.C,1)
        
        hs = halfspace(S.C(i,:),S.d(i));
        
        % intersect set with each constrained that is violated
        if ~contains_(hs,Z)
            res = and_(res,hs,'exact');
        end
    end
    
    
elseif isa(S, 'levelSet')
    
    % compute interval enclosure of constrained polynomial zonotope
    I = interval(cPZ,'interval');

    % enclose nonlinear constraint of the level set with a Taylor model
    tay = taylm(I);
    T = S.funHan(tay);

    ind = zeros(dim(I),1);
    names1 = T.names_of_var;

    for i = 1:length(ind)
        name = tay(i,1).names_of_var;
        for j = 1:length(ind)
            if strcmp(names1{j},name{1})
                ind(i) = j;
            end
        end
    end

    EC = T.monomials(2:end,2:end)';
    EC = EC(ind,:);
    A = T.coefficients(2:end)';
    b = -T.coefficients(1);

    % construct conPolyZono object for the level set
    E = eye(dim(I));
    rem = T.remainder;

    if ~strcmp(S.compOp,'==')
        temp = interval(T);
        rem = interval(infimum(temp),supremum(rem));
    end

    c = center(rem); r = rad(rem);
    if ~all(r == 0)
        A = [A,-diag(r)];
        b = b + c;
        EC = blkdiag(EC,eye(length(r)));
        E = [E;zeros(length(r),size(E,2))];
    end

    ls = conPolyZono(center(I),zonotope(I).G,E,A,b,EC);

    % intersect with the original conPolyZono object
    res = and_(cPZ,ls,'exact');
    
else
    throw(CORAerror('CORA:noops',cPZ,S));
end

% ------------------------------ END OF CODE ------------------------------
