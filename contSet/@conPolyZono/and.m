function res = and(obj,S)
% and - Computes the intersection of a conPolyZono object with an
%       other set representations
%
% Syntax:  
%    res = and(obj,S)
%
% Inputs:
%    obj - conPolyZono object
%    S - second set (supported objects: conPolyZono, halfspace, 
%                    conHyperlane, contSet objects)
%
% Outputs:
%    res - conPolyZonotope object
%
% Example: 
%    % intersection with halfspace
%    cPZ = conPolyZono([0;0],[2 0 2;0 2 2],[1 0 3;0 1 1]);
%    hs  = halfspace([1 1],0);
%
%    res1 = cPZ & hs;
%
%    figure; hold on;
%    xlim([-5,5]); ylim([-5,5]);
%    plot(hs,[1,2],'r');
%    plot(res1,[1,2],'b','EdgeColor','none','Filled',true,'Splits',20);
%    plot(cPZ,[1,2],'g','LineWidth',2);
%
%    % intersection with constrained hyperplane
%    cPZ = conPolyZono([0;0],[2 0 2;0 2 2],[1 0 3;0 1 1]);
%    ch  = conHyperplane([1 1],0,[-1 0],1);
%
%    res2 = cPZ & ch;
%
%    figure; hold on;
%    xlim([-5,5]); ylim([-5,5]);
%    plot(ch,[1,2],'r');
%    plot(res2,[1,2],'b','EdgeColor','none','Filled',true,'Splits',15);
%    plot(cPZ,[1,2],'g','LineWidth',2);
%
%    % intersection of two constrained polynomial zonotopes
%    cPZ1 = conPolyZono([0;0],[2 0 2;0 2 2],[1 0 3;0 1 1]);
%    cPZ2 = conPolyZono([-1;1],[2 0 2;0 4 -4],[1 0 2;0 1 1]);
%
%    res3 = cPZ1 & cPZ2;
%
%    figure; hold on;
%    plot(res3,[1,2],'Filled',true,'FaceColor',[0 .5 0],'EdgeColor', ...
%         'none','Splits',20);
%    plot(cPZ1,[1,2],'r','LineWidth',1.5);
%    plot(cPZ2,[1,2],'b','LineWidth',1.5);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/and, interval/and

% Author:  Niklas Kochdumper
% Written: 07-November-2018
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % determine which set is the conPolyZono object
    if ~isa(obj,'conPolyZono')
        temp = obj;
        obj = S;
        S = temp;
    end
    
    % different cases for different set represnetations
    if isa(S,'conPolyZono') || isa(S,'mptPolytope') || ...
       isa(S,'interval') || isa(S,'zonotope') || ...
       isa(S,'zonoBundle') || isa(S,'conZonotope') || ...
       isa(S,'ellipsoid') || isa(S,'capsule') || ...
       isa(S,'polyZonotope') || isa(S,'taylm')

        % convert to constrained polynomial zonotope
        S = conPolyZono(S);
        
        % compute constraint part of the resulting set
        obj.b = [obj.b; S.b; S.c - obj.c];
        A = [obj.G  -S.G];
        expMat_ = blkdiag(obj.expMat,S.expMat);
        
        if isempty(obj.A) && isempty(S.A)         
            obj.A = A;
            obj.expMat_ = expMat_;
        elseif isempty(obj.A)
            obj.A = blkdiag(S.A,A);
            O = zeros(size(obj.expMat,1),size(S.expMat_,2));
            obj.expMat_ = [[O;S.expMat_],expMat_];
        elseif isempty(S.A)
            obj.A = blkdiag(obj.A,A);
            O = zeros(size(S.expMat,1),size(obj.expMat_,2));
            obj.expMat_ = [[obj.expMat_;O],expMat_];
        else           
            obj.A = blkdiag(obj.A,S.A,A);
            O1 = zeros(size(S.expMat,1),size(obj.expMat_,2));
            O2 = zeros(size(obj.expMat,1),size(S.expMat_,2)); 
            obj.expMat_ = [[obj.expMat_;O1],[O2;S.expMat_],expMat_];
        end
        
        % compute state part of the resulting set
        O = zeros(size(S.expMat,1),size(obj.expMat,2));
        obj.expMat = [obj.expMat;O];
        
        % compute independent part of the resulting set
        if ~isempty(obj.Grest)
            if ~isempty(S.Grest)
               n = dim(obj); cen = zeros(n,1);
               zono1 = zonotope(cen,obj.Grest);
               zono2 = zonotope(cen,S.Grest);
               zono = enclose(zono1,zono2);
               obj.Grest = generators(zono);
            end
        else
            obj.Grest = S.Grest;
        end
        
        % update identifier vector
        obj.id = [obj.id; (max(obj.id)+1:max(obj.id)+length(S.id))'];
        
        % remove redundant monomials
        res = compact(obj);

        
    elseif isa(S, 'halfspace')

        % compute lower bound
        l = supportFunc(obj,S.c,'lower','interval');
        
        % add additional constraints
        A = [S.c' * obj.G, -0.5*(S.d-l)];
        b = 0.5*(S.d+l) - S.c'*obj.c;
        expMat_ = blkdiag(obj.expMat,1);
        
        obj.expMat = [obj.expMat; zeros(1,size(obj.expMat,2))];
        obj.id = [obj.id;max(obj.id) + 1];
        
        if isempty(obj.A)
           obj.A = A; obj.b = b; obj.expMat_ = expMat_;
        else
           obj.A = blkdiag(obj.A, A);
           obj.b = [obj.b; b];
           obj.expMat_ = [[obj.expMat_;zeros(1,size(obj.expMat_,2))], ...
                           expMat_];
        end
        
        % remove redundant monomials
        res = compact(obj);

        
    elseif isa(S, 'conHyperplane')

        % calculate intersection between conPolyZono object and hyperplane
        % defined as c * x = d
        obj.b = [obj.b; S.h.d - S.h.c' * obj.c];
        obj.expMat_ = [obj.expMat_ obj.expMat];
        
        if isempty(obj.A)
           obj.A = S.h.c' * obj.G; 
        else
           obj.A = blkdiag(obj.A, S.h.c' * obj.G); 
        end
        
        res = obj;

        % check if the zonotope over-approximation violates the constraints
        % (fast test to see if the computational expensive instersection 
        % with the constraits has to be performed)
        zono = zonotope(res);

        for i = 1:size(S.C,1)
            
            hs = halfspace(S.C(i,:),S.d(i));
            
            % intersect set with each constrained that is violated
            if ~in(hs,zono)
                res = res & hs;
            end
        end
        
        
    elseif isa(S, 'levelSet')
        
        % compute interval enclosure of constrained polynomial zonotope
        int = interval(obj,'interval');

        % enclose nonlinear constraint of the level set with a Taylor model
        tay = taylm(int);
        T = S.funHan(tay);

        ind = zeros(dim(int),1);
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
        expMat = eye(dim(int));
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

        ls = conPolyZono(center(int),generators(zonotope(int)),expMat, ...
                         A,b,expMat_);

        % intersect with the original conPolyZono object
        res = obj & ls;
        
    else
        error(noops(obj,S));
    end
end

%------------- END OF CODE --------------