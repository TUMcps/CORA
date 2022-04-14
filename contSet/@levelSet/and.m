function res = and(obj,S)
% and - Computes the intersection of between a levelSet and a set S with
%       the method described in [1]
%
% Syntax:  
%    res = and(obj,S)
%
% Inputs:
%    obj - level set (class: levelSet)
%    S - contSet object
%
% Outputs:
%    res - over-approximation of the intersection
%
% Example: 
%    syms x y
%    eq = x^2 + y^2 - 4;
%    ls = levelSet(eq,[x;y],'==');
%    
%    pZ = polyZonotope([2;2],[0.5 0 0.5;0 0.5 0.5],[],eye(3));
%
%    res = ls & pZ;
%
%    figure
%    hold on
%    plot(pZ,[1,2],'g','Filled',true,'EdgeColor','none');
%    xlim([-4,4]);
%    ylim([-4,4]);
%    plot(ls,[1,2],'b');   
%    plot(res,[1,2],'r','Filled',true,'EdgeColor','r','LineWidth',3);
%
% References:
%   [1] N. Kochdumper et al. "Reachability Analysis for Hybrid Systems with 
%       Nonlinear Guard Sets", HSCC 2020
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conHyperplane/and, halfspace/and

% Author: Niklas Kochdumper
% Written: 22-July-2019
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % check input arguments
    if isa(S,'conPolyZono')
        res = S & obj; return;
    end

    if ~strcmp(obj.compOp,'==')
       error('Not implemented yet!'); 
    end

    % compute interval over-approximation
    int = interval(S);
    
    % tighten interval to enclose the level set
    int = tightenDomain(obj,int);
    
    % compute polynomial zonotope with unsolvable method
    [res,err,var] = polyZonotopeUnsolvable(obj,int);
    
    % check if equation is solvable for one variable
    if obj.solveable
        
        % select variable for taylor expansion
        [var_,eq] = selectVariable(obj,int);
        
        if length(eq) == 1
            
            % compute polynomial zonotope with the solvable method
            [res_,err_] = polyZonotopeSolveable(eq{1},var_,int);
            
            % select the better over-approximation
            if err_ < err
                res = res_;
                err = err_;
                var = var_;
            end
        end
    end
    
    % use interval enclosure if it is smaller
    if err > 2*rad(int(var))
       res = polyZonotope(int); 
    end
    
    % convert back to original set representation
    if ~isa(S,'polyZonotope')
       contSet = class(S);
       eval(['res = ',contSet,'(res)';]);
    end
end


% Auxiliary Functions -----------------------------------------------------

function [res,err,var] = polyZonotopeUnsolvable(obj,int)
% compute over-approximating polynomial zonotope for the case where the
% nonlinear equation of the level set is not solvable for one variable 
% (see Sec. 4.3 in [1])
    
    % compute matrices of taylor expansion
    m = center(int);
    r = rad(int);
    n = length(m);
    int_ = int - m;

    eq = obj.funHan(m);
    grad = obj.der.grad(m);
    hess = obj.der.hess(m);

    third = cell(length(m),1);
    for k = 1:length(third)
       han = obj.der.third{k};
       third{k} = han(int);
    end
    
    % select variable for which the expansion is solved
    [~,var] = max(abs(grad));
    
    % compute lagrange remainder
    rem = interval(0,0);
    
    for i = 1:length(third)
        rem = rem + 1/6 * int_(i) * quadEval(third{i},int_);
    end
    
    % compute polynomial zonotope (center and remainder)
    c_ = -eq - center(rem) + grad(var)*m(var);
    Grest_ = rad(rem);
    
    % compute polynomial zonotope (quadratic term)
    G_ = zeros(1,n/2 + n^2/2);
    expMat_ = zeros(n,size(G_,2));
    counter = 1;
    
    for i = 1:n
        
        G_(counter) = r(i)^2 * hess(i,i);
        expMat_(i,counter) = 2;
        counter = counter + 1;
        
        for j = i+1:n
            G_(counter) = r(i) * r(j) * (hess(i,j) + hess(j,i));
            expMat_(i,counter) = 1;
            expMat_(j,counter) = 1;
            counter = counter + 1;
        end
    end
    
    % compute polynomial zonotope (linear term)
    temp = grad' * diag(r);
    temp(:,var) = [];
    G_ = -[0.5*G_,temp];
    
    temp = eye(n);
    temp(:,var) = [];
    expMat_ = [expMat_,temp];
    
    % assemble resulting polynomial zonotope
    c = zeros(n,1);
    G = zeros(n,size(G_,2) + n -1);
    Grest = zeros(n,size(Grest_,2));
    
    ind = setdiff(1:n,var);
    
    c(var) = c_./grad(var);
    c(ind) = m(ind);
    
    G(var,1:size(G_,2)) = G_./grad(var);
    G(ind,size(G_,2)+1:end) = diag(r(ind));
    
    temp = eye(n);
    temp(:,var) = [];
    expMat = [expMat_, temp];
    
    Grest(var,:) = Grest_./grad(var);
    
    res = polyZonotope(c,G,Grest,expMat);
    
    % approximation error added as uncertainty
    err = 2*rad(rem + 0.5*quadEval(hess,int-m));
end

function [res,err] = polyZonotopeSolveable(eq,var,int)
% compute over-approximating polynomial zonotope for the case where the
% nonlinear equation of the level set is solvable for one variable 
% (see Sec. 4.3 in [1])

    % interval without selected variable
    infi = infimum(int);
    sup = supremum(int);
    
    infi(var) = [];
    sup(var) = [];
    
    int_ = interval(infi,sup);
    m_ = center(int_);
    r_ = rad(int_);

    % taylor series expansion point == center of interval
    m = center(int);
    n = length(m);
    
    % compute matrices of taylor expansion
    f = eq.eq(m);
    grad = eq.grad(m);
    hess = eq.hess(m);
    third = cell(length(eq.third),1);
    for i = 1:length(third)
       funHan = eq.third{i};
       third{i} = funHan(int);
    end
    
    % compute lagrange remainder
    intTemp = int_ - m_;
    rem = interval(0,0);
    
    for i = 1:length(third)
        rem = rem + 1/6 * intTemp(i) * quadEval(third{i},intTemp);
    end

    % compute polynomial zonotope (center and remainder)
    c_ = f + center(rem);
    Grest_ = rad(rem);
    
    % compute polynomial zonotope (quadratic term)
    G_ = zeros(1,(n-1)/2 + (n-1)^2/2);
    expMat_ = zeros(n-1,size(G_,2));
    counter = 1;
    
    for i = 1:n-1
        
        G_(counter) = r_(i)^2 * hess(i,i);
        expMat_(i,counter) = 2;
        counter = counter + 1;
        
        for j = i+1:n-2
            G_(counter) = r_(i) * r_(j) * (hess(i,j) + hess(j,i));
            expMat_(i,counter) = 1;
            expMat_(j,counter) = 1;
            counter = counter + 1;
        end
    end
    
    % compute polynomial zonotope (linear term)
    G_ = [0.5*G_,grad' * diag(r_)];
    expMat_ = [expMat_,eye(n-1)];
    
    
    % assemble polynomial zonotope
    if var == 1
       c = [c_;m_];
       Grest = [Grest_;zeros(n-1,1)];
       G = blkdiag(G_,diag(r_));
    elseif var == n
       c = [m_;c_];
       Grest = [zeros(n-1,1);Grest_];
       G = [[zeros(n-1,size(G_,2));G_],[diag(r_);zeros(1,n-1)]];
    else
       c = [m(1:var-1);c_;m(var+1:end)];
       Grest = [zeros(var-1,1);Grest_;zeros(n-var,1)];
       temp = [r_(1:var-1);0;r_(var:end)];
       temp = diag(temp);
       temp(:,var) = [];
       G = [[zeros(var-1,size(G_,2));G_;zeros(n-var,size(G_,2))],temp];
    end
    
    expMat = [expMat_,eye(n-1)];
    
    res = polyZonotope(c,G,Grest,expMat);
    
    % approximation error
    err = 2*rad(rem);
    
end

function [var,eq] = selectVariable(obj,int)
% select variable for taylor expansion

    % check if there is one equation with a unique solution
    for i = 1:length(obj.solved)
        if obj.solved{i}.contained && obj.solved{i}.solveable && length(obj.solved{i}.eq) == 1
           var = i;
           eq = {obj.solved{i}.funHan{1}};
           return;
        end
    end

    % select the best equation from the ones that have two solutions
    var = [];
    eq = [];
    
    infi = infimum(int);
    sup = supremum(int);

    for i = 1:length(obj.solved)
        
       if obj.solved{i}.contained && obj.solved{i}.solveable
           
          eqTemp = {};
          
          for j = 1:length(obj.solved{i}.eq)
              
              funHan = obj.solved{i}.funHan{j}.eq;
              
              try
                  intTemp = funHan(int);
                  if ~(supremum(intTemp) < infi(i) || ...
                       infimum(intTemp) > sup(i))

                        eqTemp{end+1} = obj.solved{i}.funHan{j};
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

function res = quadEval(Q,int)
% tight evaluation of a quadratic term using interval arithmetic

    res = interval(0,0);
        
    for k = 1:length(int)
        temp = interval(0,0);
        for l = k+1:length(int)
            temp = temp + (Q(k,l) + Q(l,k)) * int(l);
        end
        res = res + Q(k,k) * int(k)^2 + temp * int(k);
    end
end

%------------- END OF CODE --------------