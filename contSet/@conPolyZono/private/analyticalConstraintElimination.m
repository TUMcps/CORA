function listFin = analyticalConstraintElimination(obj)  
% analyticalConstraintElimination - Eliminate some special cases of
%                                   constraints that can be solved
%                                   analytically
%
% Syntax:  
%    listFin = analyticalConstraintElimination(obj)  
%
% Inputs:
%    obj - conPolyZonotope object
%
% Outputs:
%    listFin - list of sets whose union is identical to the original set
%
% Example: 
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: plot, polygon

% Author:       Niklas Kochdumper
% Written:      30-November-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    list{1} = obj;
    listFin = {};
    
    % loop until no more constraints can be eliminated
    while ~isempty(list)
        
        listNew = {}; 
        
        for i = 1:length(list)
           
            objTemp = list{i};
            val = [];
            indices = [];
            
            % constraints with only one single entry
            ind = find(sum(objTemp.A ~= 0,2) == 1);
            
            for j = 1:length(ind)
               
                [val_,ind_,empty] = elimSingleEntryConstraint( ...
                                    objTemp.A(ind(j),:),objTemp.b(ind(j)),objTemp.expMat_);
                
                if empty
                    continue;
                elseif ~isempty(val_)
                
                    [val,indices] = updateValues(val,indices,val_,ind_);

                    if isempty(val)
                       continue; 
                    end
                end
            end
           
            % quadratic constraints
            for j = 1:size(objTemp.A,1)
               ind = find(objTemp.A(j,:) ~= 0);
               if max(sum(objTemp.expMat_(:,ind),1)) == 2
                   [val_,ind_,empty] = elimQuadConstraint(objTemp.A(j,ind), ...
                                            objTemp.b(j),objTemp.expMat_(:,ind));
                                        
                   if empty
                      continue;
                   elseif ~isempty(val_)
                       
                      [val,indices] = updateValues(val,indices,val_,ind_);

                      if isempty(val)
                         continue; 
                      end
                   end
               end
            end
            
            % special case of constraint
            for j = 1:size(objTemp.A,1)
                if objTemp.b(j) == 0
                    ind = find(objTemp.A(j,:) ~= 0);
                    if all(objTemp.expMat_(:,ind) == 0 | objTemp.expMat_(:,ind) == 2)
                        [val_,ind_] = elimSpecialCase(objTemp.A(ind),objTemp.expMat_(:,ind));
                        
                        if ~isempty(val_)
                            [val,indices] = updateValues(val,indices,val_,ind_);

                             if isempty(val)
                                continue; 
                             end
                        end
                    end
                end
            end
            
            % substitute the values for the factors into the set
            if ~isempty(indices)
                
                % insert all possible value combinations into the object
                listTemp = cell(size(val,2),1);
                
                for j = 1:length(listTemp)
                    listTemp{j} = subsFactorValue(objTemp,indices,val(:,j));
                end
                
                 % update list
                 listNew = [listNew;listTemp];
                
            else
                listFin{end+1} = objTemp;
            end
        end
        
        % update list
        list = listNew;
        
    end
end




% Auxiliary Functions -----------------------------------------------------

function [val_,ind_,empty] = elimSingleEntryConstraint(A,b,expMat)

    empty = 0; val_ = [];

    index = find(A);
    ind_ = find(expMat(:,index));

    % fix a single factor
    if length(ind_) == 1

        temp = nthroot(b/A(index),expMat(ind_,index));

        if temp > 1     % value located outside the domain
           empty = 1;
           return;
        else
           if mod(expMat(ind_,index),2) == 0
               val_ = [-temp,temp];
           else
               val_ = temp;
           end
        end

    % fix two factors
    elseif length(ind_) == 2

        if expMat(ind_(1),index) == 1 && expMat(ind_(2),index) == 1

           if b/A(index) == 1
               val_ = [1 -1;1 -1];
           elseif b/A(index) == -1
               val_ = [1 -1;-1 1];
           elseif b/A(index) > 1 || bA(index) < -1
               empty = 1;
               return;
           end
        end
    end 
end
    
function [val,ind,empty] = elimQuadConstraint(A,b,expMat)

    empty = 0;
    val = [];
    
    % eliminate zero rows from the exponent matrix
    ind = find(sum(expMat,2) > 0);
    expMat = expMat(ind,:);

    % check for positive or negative definitness
    H = hessian(A,expMat); 
    [~,p] = chol(H);
    
    if p == 0   % Hessian matrix positive definite

        % compute point where the function reaches its maximum
        l = linPart(A,expMat);
        x = linsolve(H,l);
        
        % check if maximum is located inside the domain
        if all(x > -1) && all(x < 1)
            
            % compute function value at the maximum
            v = l'*x + x' * H * x;
            
            % check if the constraint can be eliminated
            if b > v
               ind = [];
               empty = 1;
            elseif b == v
               val = x;
            end
        end
        
    else
        [~,p] = chol(-H);
        
        if p == 0   % Hessian matrix negative definite
            
            % compute point where the function reaches its minimum
            l = linPart(A,expMat);
            x = linsolve(H,l);

            % check if minimum is located inside the domain
            if all(x > -1) && all(x < 1)

                % compute function value at the minimum
                v = l'*x + x' * H * x;

                % check if the constraint can be eliminated
                if b < v
                   ind = [];
                   empty = 1;
                elseif b == v
                   val = x;
                end
            end
        end
    end
end


function [val,ind] = elimSpecialCase(A,expMat)

    val = [];

    % substitue all square terms by new variables
    expMat_ = double(expMat > 0);
    
    % compute the gradients at [0;0; ... ;0]
    grad = zeros(size(expMat,1),1);
    z = zeros(length(grad),1);
    
    for i = 1:length(grad)
       indTemp = find(expMat_(i,:) > 0);
       Etemp = expMat_(:,indTemp);
       Etemp(i,:) = zeros(1,size(Etemp,2));
       Atemp = A(indTemp);
       
       for j = 1:size(Etemp,2)
          grad(i) = grad(i) + Atemp(j) * prod(z.^Etemp(:,j)); 
       end
    end
    
    % determine variables with gradient greater than 0
    ind = find(grad > 0);
    
    if ~isempty(ind)
        % evaluate gradient on the whole domain using interval arithmetic
        temp = ones(size(expMat,1),1);
        dom = interval(0*temp,temp);
        gradInt = interval(zeros(length(ind),1));

        for i = 1:length(ind)
           indTemp = find(expMat_(i,:) > 0);
           Etemp = expMat_(:,indTemp);
           Etemp(i,:) = zeros(1,size(Etemp,2));
           Atemp = A(indTemp);        

           for j = 1:size(Etemp,2)
              product = dom(1).^Etemp(1,j);
              for k = 2:length(dom)
                 product = product * dom(k).^Etemp(k,j); 
              end

              gradInt(i) = gradInt(i) + Atemp(j) * product;
           end
        end

        % check if minimum of the function is reached at [0;0; ... ;0]
        if all(infimum(gradInt) > 0)
           val = zeros(length(ind),1);
        end
    end
end


function H = hessian(A,expMat)

    H = zeros(size(expMat,1));
    for i = 1:size(expMat,2)
       ind = find(expMat(:,i) > 0);
       if length(ind) == 2
          H(ind(1),ind(2)) = A(i);
          H(ind(2),ind(1)) = A(i);
       else
          H(ind(1),ind(1)) = A(i); 
       end
    end
end

function b = linPart(A,expMat)

    % find monomials with polynomial degree 1
    ind = find(sum(expMat,1) == 1);
    
    % determine slope for each factor
    b = zeros(size(expMat,1),1);
    
    for i = 1:length(b)
        indTemp = find(expMat(i,ind) > 0);
        if ~isempty(indTemp)
           b(i) = sum(A(:,ind(indTemp)),2); 
        end
    end
end

function [val,indices] = updateValues(val,indices,val_,ind_)

    % no stored values exist
    if isempty(indices)
        val = val_;
        indices = ind_;
        
    % merge with the stored values
    else
        
        [int,ind1,ind2] = intersect(indices,ind_);
        
        % no common factor between stored and new values
        if isempty(int)
           temp = ones(1,size(val,2));
           valOld = val;
           val = [];
           for i = 1:size(val_,2)
                val = [val, [val_(:,i)*temp ; valOld]];
           end
           indices = [ind_;indices];
           
        % merge the values for the common factors
        else
            
           % merge indices 
           ind2rem = setdiff(1:length(ind_),ind2);           
           indicesNew = [indices;ind_(ind2rem)];
            
           % merge values
           valNew = [];
           
           for i = 1:size(val_,2)
               index = find(all(bsxfun(@eq, val_(ind2,i), val(ind1,:)), 1));  
               for j = 1:length(index)
                    valNew = [valNew, [val(:,index(j));val_(ind2rem,i)]];
               end
           end
           
           % assign output arguments
           val = valNew;
           indices = indicesNew;
           
        end
    end
end

function res = subsFactorValue(obj,ind,val)
% substitute factors with a numerical values

    % extract object properties
    expMat = obj.expMat;
    G = obj.G;
    expMat_ = obj.expMat_;
    A = obj.A;

    % substite value into generators
    temp = val.^obj.expMat(ind,:);
    
    for i = 1:size(temp,1)
        G = G * diag(temp(i,:));
    end
    
    expMat(ind,:) = [];
    
    % substite value into constraints
    temp = val.^obj.expMat_(ind,:);
    
    for i = 1:size(temp,1)
        A = A * diag(temp(i,:));
    end
    
    expMat_(ind,:) = [];
    
    % remove factor identifier
    id = obj.id;
    id(ind) = [];
    
    % construct resulting conPolyZonotope object
    res = conPolyZono(obj.c,G,expMat,A,obj.b,expMat_,obj.Grest,id);
    
    % remove redundant monomials
    res = compact(res);
    
    % remove trivial constriants
    res = removeTrivialConstraints(res);

end

function res = removeTrivialConstraints(obj)
% remove the trivial constraint 0 = 0*a1 + 0*a2 + ...

    temp = sum(abs(obj.A),2);
    
    A = obj.A;
    b = obj.b;
    ind = [];
    
    for i = 1:length(temp)
       if temp(i) == 0 && obj.b(i) == 0
           ind = [ind;i];
       end
    end

    if ~isempty(ind)
        A(ind,:) = [];
        b(ind) = [];
        res = conPolyZono(obj.c,obj.G,obj.expMat,A,b,obj.expMat_, ...
                          obj.Grest,obj.id);
    else
        res = obj; 
    end
end
    
%------------- END OF CODE --------------