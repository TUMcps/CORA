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
%    plot(res1,[1,2],'FaceColor','b','Splits',20);
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
%    plot(res2,[1,2],'FaceColor','b','Splits',15);
%    plot(cPZ,[1,2],'g','LineWidth',2);
%
%    % intersection of two constrained polynomial zonotopes
%    cPZ1 = conPolyZono([0;0],[2 0 2;0 2 2],[1 0 3;0 1 1]);
%    cPZ2 = conPolyZono([-1;1],[2 0 2;0 4 -4],[1 0 2;0 1 1]);
%
%    res3 = cPZ1 & cPZ2;
%
%    figure; hold on;
%    plot(res3,[1,2],'FaceColor',[0 .5 0],'Splits',20);
%    plot(cPZ1,[1,2],'r','LineWidth',1.5);
%    plot(cPZ2,[1,2],'b','LineWidth',1.5);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/and_, interval/and_

% Author:       Niklas Kochdumper
% Written:      07-November-2018
% Last update:  ---
% Last revision:27-March-2023 (MW, rename and_)

%------------- BEGIN CODE --------------
           
% different cases for different set representations
if isa(S,'conPolyZono') || isa(S,'mptPolytope') || ...
   isa(S,'interval') || isa(S,'zonotope') || ...
   isa(S,'zonoBundle') || isa(S,'conZonotope') || ...
   isa(S,'ellipsoid') || isa(S,'capsule') || ...
   isa(S,'polyZonotope') || isa(S,'taylm')

    % convert to constrained polynomial zonotope
    S = conPolyZono(S);
    
    % compute constraint part of the resulting set
    cPZ.b = [cPZ.b; S.b; S.c - cPZ.c];
    A = [cPZ.G  -S.G];
    expMat_ = blkdiag(cPZ.expMat,S.expMat);
    
    if isempty(cPZ.A) && isempty(S.A)         
        cPZ.A = A;
        cPZ.expMat_ = expMat_;
    elseif isempty(cPZ.A)
        cPZ.A = blkdiag(S.A,A);
        O = zeros(size(cPZ.expMat,1),size(S.expMat_,2));
        cPZ.expMat_ = [[O;S.expMat_],expMat_];
    elseif isempty(S.A)
        cPZ.A = blkdiag(cPZ.A,A);
        O = zeros(size(S.expMat,1),size(cPZ.expMat_,2));
        cPZ.expMat_ = [[cPZ.expMat_;O],expMat_];
    else           
        cPZ.A = blkdiag(cPZ.A,S.A,A);
        O1 = zeros(size(S.expMat,1),size(cPZ.expMat_,2));
        O2 = zeros(size(cPZ.expMat,1),size(S.expMat_,2)); 
        cPZ.expMat_ = [[cPZ.expMat_;O1],[O2;S.expMat_],expMat_];
    end
    
    % compute state part of the resulting set
    O = zeros(size(S.expMat,1),size(cPZ.expMat,2));
    cPZ.expMat = [cPZ.expMat;O];
    
    % compute independent part of the resulting set
    if ~isempty(cPZ.Grest)
        if ~isempty(S.Grest)
           n = dim(cPZ); cen = zeros(n,1);
           Z1 = zonotope(cen,cPZ.Grest);
           Z2 = zonotope(cen,S.Grest);
           Z = enclose(Z1,Z2);
           cPZ.Grest = generators(Z);
        end
    else
        cPZ.Grest = S.Grest;
    end
    
    % update identifier vector
    cPZ.id = [cPZ.id; (max(cPZ.id)+1:max(cPZ.id)+length(S.id))'];
    
    % remove redundant monomials
    res = compact(cPZ);

    
elseif isa(S,'halfspace')

    % compute lower bound
    l = supportFunc_(cPZ,S.c,'lower','interval',[]);
    
    % add additional constraints
    A = [S.c' * cPZ.G, -0.5*(S.d-l)];
    b = 0.5*(S.d+l) - S.c'*cPZ.c;
    expMat_ = blkdiag(cPZ.expMat,1);
    
    cPZ.expMat = [cPZ.expMat; zeros(1,size(cPZ.expMat,2))];
    cPZ.id = [cPZ.id;max(cPZ.id) + 1];
    
    if isempty(cPZ.A)
       cPZ.A = A; cPZ.b = b; cPZ.expMat_ = expMat_;
    else
       cPZ.A = blkdiag(cPZ.A, A);
       cPZ.b = [cPZ.b; b];
       cPZ.expMat_ = [[cPZ.expMat_;zeros(1,size(cPZ.expMat_,2))], ...
                       expMat_];
    end
    
    % remove redundant monomials
    res = compact(cPZ);

    
elseif isa(S, 'conHyperplane')

    % calculate intersection between conPolyZono object and hyperplane
    % defined as c * x = d
    cPZ.b = [cPZ.b; S.h.d - S.h.c' * cPZ.c];
    cPZ.expMat_ = [cPZ.expMat_ cPZ.expMat];
    
    if isempty(cPZ.A)
       cPZ.A = S.h.c' * cPZ.G; 
    else
       cPZ.A = blkdiag(cPZ.A, S.h.c' * cPZ.G); 
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

    expMat_ = T.monomials(2:end,2:end)';
    expMat_ = expMat_(ind,:);
    A = T.coefficients(2:end)';
    b = -T.coefficients(1);

    % construct conPolyZono object for the level set
    expMat = eye(dim(I));
    rem = T.remainder;

    if ~strcmp(S.compOp,'==')
        temp = interval(T);
        rem = interval(infimum(temp),supremum(rem));
    end

    c = center(rem); r = rad(rem);
    if ~all(r == 0)
        A = [A,-diag(r)];
        b = b + c;
        expMat_ = blkdiag(expMat_,eye(length(r)));
        expMat = [expMat;zeros(length(r),size(expMat,2))];
    end

    ls = conPolyZono(center(I),generators(zonotope(I)),expMat, ...
                     A,b,expMat_);

    % intersect with the original conPolyZono object
    res = and_(cPZ,ls,'exact');
    
else
    throw(CORAerror('CORA:noops',cPZ,S));
end

%------------- END OF CODE --------------