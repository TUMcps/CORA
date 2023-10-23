function res = and_(ls,S,varargin)
% and_ - Computes the intersection of between a level set and a set S with
%    the method described in [1]
%
% Syntax:
%    res = and_(ls,S)
%
% Inputs:
%    ls - levelSet object
%    S - contSet object
%
% Outputs:
%    res - over-approximation of the intersection
%
% Example: 
%    syms x y
%    eq = x^2 + y^2 - 4;
%    ls = levelSet(eq,[x;y],'==');
%    pZ = polyZonotope([2;2],[0.5 0 0.5;0 0.5 0.5],[],eye(3));
%
%    res = ls & pZ;
%
%    figure; hold on; xlim([-4,4]); ylim([-4,4]);
%    plot(pZ,[1,2],'FaceColor','g');
%    plot(ls,[1,2],'b');   
%    plot(res,[1,2],'FaceColor','r','LineWidth',3);
%
% References:
%   [1] N. Kochdumper et al. "Reachability Analysis for Hybrid Systems with 
%       Nonlinear Guard Sets", HSCC 2020
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and, conHyperplane/and_, halfspace/and_

% Authors:       Niklas Kochdumper
% Written:       22-July-2019
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename and_)

% ------------------------------ BEGIN CODE -------------------------------

    % ensure order
    [ls,S] = findClassArg(ls,S,'levelSet');

    % intersection with fullspace or emptySet
    if isa(S,'fullspace') || isa(S,'emptySet')
        res = and_(S,ls); return
    end

    % check input arguments
    if isa(S,'conPolyZono')
        res = and_(S,ls,'exact'); return
    end

    % convert polytope to level set
    if isa(S,'polytope')
        S = levelSet(S);
    end

    % special case: intersection with a point
    if isa(S,'polyZonotope') && aux_isavector(S)
        % check if point fulfills level set equation
        if ls.funHan(S.c) < 0 || withinTol(ls.funHan(S.c),0)
            % init as zonotope for upcoming computations
            res = zonotope(S.c);
        else
            res = emptySet(dim(S));
        end
        return
    end
    
    % check if S is a level Set with same compOp => concatenate
    if isa(S,'levelSet')
        % use vars from ls (should be irrelevant which ones are used)
        vars = ls.vars;
        newEqs = [ls.eq;subs(S.eq,S.vars,vars)];
        newCompOp = aux_uniteCompOp(ls.compOp,S.compOp);
        res = levelSet(newEqs,vars,newCompOp);
        res = compact_(res,'all',eps);
        return;
    end

    % split up ls if it has multiple equations
    if iscell(ls.compOp) && not(isscalar(ls.compOp))
        res = fullspace(length(ls.vars));
        for i = 1:length(ls.compOp)
            op = ls.compOp{i};
            eq = ls.eq(i);
            lsComp = levelSet(eq, ls.vars, op);
            comp = lsComp & S;
            res = res & comp;
        end
        return;
    end

    % compute coarse outer-approximation
    if any(strcmp(ls.compOp,{'<=','<'}))
        % caution: multiple equations not supported...

        % outer-approximation second set by interval for range bounding
        I = interval(S);

        % perform range bounding using Taylor models for complement of
        % level set and second set
        ls_ = not(ls);
        boundedVals = interval(taylm(symfun(ls_.eq, ls_.vars),I));

        % if the entire range is contained in the complement, the
        % intersection is empty
        % (TODO: slightly different for '<' / '<=')
        if contains_(interval(-Inf,0),boundedVals)
            res = emptySet(dim(S));
            return
        end

        % S is an outer-approximation of S & ls, use this for now
        % TODO: implement contractors for tighter outer-approximation
        res = zonotope(S);
        return

    end
    
    if ~strcmp(ls.compOp,'==')
        throw(CORAerror('CORA:noops',ls,S));
    end

    % compute interval over-approximation
    I = interval(S);
    
    % tighten interval to enclose the level set
    I = tightenDomain(ls,I);
    
    % compute polynomial zonotope with unsolvable method
    [res,err,var] = aux_polyZonotopeUnsolvable(ls,I);
    
    % check if equation is solvable for one variable
    if ls.solvable
        
        % select variable for taylor expansion
        [var_,eq] = aux_selectVariable(ls,I);
        
        if length(eq) == 1
            
            % compute polynomial zonotope with the solvable method
            [res_,err_] = aux_polyZonotopeSolvable(eq{1},var_,I);
            
            % select the better over-approximation
            if err_ < err
                res = res_;
                err = err_;
                var = var_;
            end
        end
    end
    
    % use interval enclosure if it is smaller
    if err > 2*rad(I(var))
        res = polyZonotope(I); 
    end
    
    % convert back to original set representation
    if ~isa(S,'polyZonotope')
        contSet = class(S);
        eval(['res = ',contSet,'(res);']);
    end
end


% Auxiliary functions -----------------------------------------------------

function [res,err,var] = aux_polyZonotopeUnsolvable(ls,I)
% compute over-approximating polynomial zonotope for the case where the
% nonlinear equation of the level set is not solvable for one variable 
% (see Sec. 4.3 in [1])
    
    % compute matrices of taylor expansion
    m = center(I);
    r = rad(I);
    n = length(m);
    int_ = I - m;

    eq = ls.funHan(m);
    grad = ls.der.grad(m);
    hess = ls.der.hess(m);

    third = cell(length(m),1);
    for k = 1:length(third)
       han = ls.der.third{k};
       third{k} = han(I);
    end
    
    % select variable for which the expansion is solved
    [~,var] = max(abs(grad));
    
    % compute lagrange remainder
    rem = interval(0,0);
    
    for i = 1:length(third)
        rem = rem + 1/6 * int_(i) * aux_quadEval(third{i},int_);
    end
    
    % compute polynomial zonotope (center and remainder)
    c_ = -eq - center(rem) + grad(var)*m(var);
    GI_ = rad(rem);
    
    % compute polynomial zonotope (quadratic term)
    G_ = zeros(1,n/2 + n^2/2);
    E_ = zeros(n,size(G_,2));
    counter = 1;
    
    for i = 1:n
        
        G_(counter) = r(i)^2 * hess(i,i);
        E_(i,counter) = 2;
        counter = counter + 1;
        
        for j = i+1:n
            G_(counter) = r(i) * r(j) * (hess(i,j) + hess(j,i));
            E_(i,counter) = 1;
            E_(j,counter) = 1;
            counter = counter + 1;
        end
    end
    
    % compute polynomial zonotope (linear term)
    temp = grad' * diag(r);
    temp(:,var) = [];
    G_ = -[0.5*G_,temp];
    
    temp = eye(n);
    temp(:,var) = [];
    E_ = [E_,temp];
    
    % assemble resulting polynomial zonotope
    c = zeros(n,1);
    G = zeros(n,size(G_,2) + n -1);
    GI = zeros(n,size(GI_,2));
    
    ind = setdiff(1:n,var);
    
    c(var) = c_./grad(var);
    c(ind) = m(ind);
    
    G(var,1:size(G_,2)) = G_./grad(var);
    G(ind,size(G_,2)+1:end) = diag(r(ind));
    
    temp = eye(n);
    temp(:,var) = [];
    E = [E_, temp];
    
    GI(var,:) = GI_./grad(var);
    
    res = polyZonotope(c,G,GI,E);
    
    % approximation error added as uncertainty
    err = 2*rad(rem + 0.5*aux_quadEval(hess,I-m));
end

function [res,err] = aux_polyZonotopeSolvable(eq,var,I)
% compute over-approximating polynomial zonotope for the case where the
% nonlinear equation of the level set is solvable for one variable 
% (see Sec. 4.3 in [1])

    % interval without selected variable
    lb = infimum(I);
    ub = supremum(I);
    
    lb(var) = [];
    ub(var) = [];
    
    I_ = interval(lb,ub);
    m_ = center(I_);
    r_ = rad(I_);

    % taylor series expansion point == center of interval
    m = center(I);
    n = length(m);
    
    % compute matrices of taylor expansion
    f = eq.eq(m);
    grad = eq.grad(m);
    hess = eq.hess(m);
    third = cell(length(eq.third),1);
    for i = 1:length(third)
       funHan = eq.third{i};
       third{i} = funHan(I);
    end
    
    % compute lagrange remainder
    intTemp = I_ - m_;
    rem = interval(0,0);
    
    for i = 1:length(third)
        rem = rem + 1/6 * intTemp(i) * aux_quadEval(third{i},intTemp);
    end

    % compute polynomial zonotope (center and remainder)
    c_ = f + center(rem);
    GI_ = rad(rem);
    
    % compute polynomial zonotope (quadratic term)
    G_ = zeros(1,(n-1)/2 + (n-1)^2/2);
    E_ = zeros(n-1,size(G_,2));
    counter = 1;
    
    for i = 1:n-1
        
        G_(counter) = r_(i)^2 * hess(i,i);
        E_(i,counter) = 2;
        counter = counter + 1;
        
        for j = i+1:n-2
            G_(counter) = r_(i) * r_(j) * (hess(i,j) + hess(j,i));
            E_(i,counter) = 1;
            E_(j,counter) = 1;
            counter = counter + 1;
        end
    end
    
    % compute polynomial zonotope (linear term)
    G_ = [0.5*G_,grad' * diag(r_)];
    E_ = [E_,eye(n-1)];
    
    % assemble polynomial zonotope
    if var == 1
       c = [c_;m_];
       GI = [GI_;zeros(n-1,1)];
       G = blkdiag(G_,diag(r_));
    elseif var == n
       c = [m_;c_];
       GI = [zeros(n-1,1);GI_];
       G = [[zeros(n-1,size(G_,2));G_],[diag(r_);zeros(1,n-1)]];
    else
       c = [m(1:var-1);c_;m(var+1:end)];
       GI = [zeros(var-1,1);GI_;zeros(n-var,1)];
       temp = [r_(1:var-1);0;r_(var:end)];
       temp = diag(temp);
       temp(:,var) = [];
       G = [[zeros(var-1,size(G_,2));G_;zeros(n-var,size(G_,2))],temp];
    end
    
    E = [E_,eye(n-1)];
    
    res = polyZonotope(c,G,GI,E);
    
    % approximation error
    err = 2*rad(rem);
    
end

function [var,eq] = aux_selectVariable(ls,I)
% select variable for taylor expansion

    % check if there is one equation with a unique solution
    for i = 1:length(ls.solved)
        if ls.solved{i}.contained && ls.solved{i}.solvable ...
                && length(ls.solved{i}.eq) == 1
           var = i;
           eq = {ls.solved{i}.funHan{1}};
           return;
        end
    end

    % select the best equation from the ones that have two solutions
    var = [];
    eq = [];
    
    lb = infimum(I);
    ub = supremum(I);

    for i = 1:length(ls.solved)
        
       if ls.solved{i}.contained && ls.solved{i}.solvable
           
          eqTemp = {};
          
          for j = 1:length(ls.solved{i}.eq)
              
              funHan = ls.solved{i}.funHan{j}.eq;
              
              try
                  intTemp = funHan(I);
                  if ~(supremum(intTemp) < lb(i) || ...
                       infimum(intTemp) > ub(i))

                        eqTemp{end+1} = ls.solved{i}.funHan{j};
                  end
              end
          end
          
          if ~isempty(eqTemp)
              if isempty(eq) || length(eq) > length(eqTemp)
                 var = i;
                 eq = eqTemp;
                 if length(eqTemp) == 1
                    break; 
                 end
              end
          end
       end
    end
end

function res = aux_quadEval(Q,I)
% tight evaluation of a quadratic term using interval arithmetic

    res = interval(0,0);
        
    for k = 1:length(I)
        temp = interval(0,0);
        for l = k+1:length(I)
            temp = temp + (Q(k,l) + Q(l,k)) * I(l);
        end
        res = res + Q(k,k) * I(k)^2 + temp * I(k);
    end
end

function compOp = aux_uniteCompOp(compOp1,compOp2)

% make all cells and vertical
if iscell(compOp1)
    if size(compOp1,2) > 1
        compOp1 = compOp1';
    end
else
    compOp1 = {compOp1};
end
if iscell(compOp2)
    if size(compOp2,2) > 1
        compOp2 = compOp2';
    end
else
    compOp2 = {compOp2};
end

% concatenate
compOp = [compOp1;compOp2];

end

function res = aux_isavector(S)
% checks if the set S represents just a vector
% TODO: integrate this function to 'representsa'
res = false;

if isa(S,'polyZonotope')
    % just a point if there are no independent generators (or all-zero) and
    % either no dependent generators (or all-zero) or the exponent matrix
    % is all-zero
    if ~isempty(S.GI) && any(any(S.GI))
        res = false; return
    elseif ~isempty(S.G) && any(any(S.G)) && any(any(S.E))
        res = false; return
    end
    res = true;
end   

end

% ------------------------------ END OF CODE ------------------------------
