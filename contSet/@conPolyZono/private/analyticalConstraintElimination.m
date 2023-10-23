function listFin = analyticalConstraintElimination(cPZ)  
% analyticalConstraintElimination - Eliminates some special cases of
%    constraints that can be solved analytically
%
% Syntax:
%    listFin = analyticalConstraintElimination(cPZ)  
%
% Inputs:
%    cPZ - conPolyZono object
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
% See also: ---

% Authors:       Niklas Kochdumper
% Written:       30-November-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    list{1} = cPZ;
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
               
                [val_,ind_,empty] = aux_elimSingleEntryConstraint(...
                    objTemp.A(ind(j),:),objTemp.b(ind(j)),objTemp.EC);
                
                if empty
                    continue;
                elseif ~isempty(val_)
                
                    [val,indices] = aux_updateValues(val,indices,val_,ind_);

                    if isempty(val)
                       continue; 
                    end
                end
            end
           
            % quadratic constraints
            for j = 1:size(objTemp.A,1)
               ind = find(objTemp.A(j,:) ~= 0);
               if max(sum(objTemp.EC(:,ind),1)) == 2
                   [val_,ind_,empty] = aux_elimQuadConstraint(objTemp.A(j,ind),...
                    objTemp.b(j),objTemp.EC(:,ind));
                                        
                   if empty
                      continue;
                   elseif ~isempty(val_)
                       
                      [val,indices] = aux_updateValues(val,indices,val_,ind_);

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
                    if all(objTemp.EC(:,ind) == 0 | objTemp.EC(:,ind) == 2)
                        [val_,ind_] = aux_elimSpecialCase(objTemp.A(ind),...
                            objTemp.EC(:,ind));
                        
                        if ~isempty(val_)
                            [val,indices] = aux_updateValues(val,indices,val_,ind_);

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
                    listTemp{j} = aux_subsFactorValue(objTemp,indices,val(:,j));
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


% Auxiliary functions -----------------------------------------------------

function [val_,ind_,empty] = aux_elimSingleEntryConstraint(A,b,E)

    empty = false; val_ = [];

    index = find(A);
    ind_ = find(E(:,index));

    % fix a single factor
    if length(ind_) == 1

        temp = nthroot(b/A(index),E(ind_,index));

        if temp > 1     % value located outside the domain
           empty = true;
           return;
        else
           if mod(E(ind_,index),2) == 0
               val_ = [-temp,temp];
           else
               val_ = temp;
           end
        end

    % fix two factors
    elseif length(ind_) == 2

        if E(ind_(1),index) == 1 && E(ind_(2),index) == 1

           if b/A(index) == 1
               val_ = [1 -1;1 -1];
           elseif b/A(index) == -1
               val_ = [1 -1;-1 1];
           elseif b/A(index) > 1 || b/A(index) < -1
               empty = true;
               return;
           end
        end
    end 
end
    
function [val,ind,empty] = aux_elimQuadConstraint(A,b,E)

    empty = false;
    val = [];
    
    % eliminate zero rows from the exponent matrix
    ind = find(sum(E,2) > 0);
    E = E(ind,:);

    % check for positive or negative definitness
    H = aux_hessian(A,E); 
    [~,p] = chol(H);
    
    if p == 0   % Hessian matrix positive definite

        % compute point where the function reaches its maximum
        l = aux_linPart(A,E);
        x = linsolve(H,l);
        
        % check if maximum is located inside the domain
        if all(x > -1) && all(x < 1)
            
            % compute function value at the maximum
            v = l'*x + x' * H * x;
            
            % check if the constraint can be eliminated
            if b > v
               ind = [];
               empty = true;
            elseif b == v
               val = x;
            end
        end
        
    else
        [~,p] = chol(-H);
        
        if p == 0   % Hessian matrix negative definite
            
            % compute point where the function reaches its minimum
            l = aux_linPart(A,E);
            x = linsolve(H,l);

            % check if minimum is located inside the domain
            if all(x > -1) && all(x < 1)

                % compute function value at the minimum
                v = l'*x + x' * H * x;

                % check if the constraint can be eliminated
                if b < v
                   ind = [];
                   empty = true;
                elseif b == v
                   val = x;
                end
            end
        end
    end
end


function [val,ind] = aux_elimSpecialCase(A,E)

    val = [];

    % substitue all square terms by new variables
    EC = double(E > 0);
    
    % compute the gradients at [0;0; ... ;0]
    grad = zeros(size(E,1),1);
    z = zeros(length(grad),1);
    
    for i = 1:length(grad)
       indTemp = find(EC(i,:) > 0);
       Etemp = EC(:,indTemp);
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
        temp = ones(size(E,1),1);
        dom = interval(0*temp,temp);
        gradInt = interval(zeros(length(ind),1));

        for i = 1:length(ind)
           indTemp = find(EC(i,:) > 0);
           Etemp = EC(:,indTemp);
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


function H = aux_hessian(A,E)

    H = zeros(size(E,1));
    for i = 1:size(E,2)
       ind = find(E(:,i) > 0);
       if length(ind) == 2
          H(ind(1),ind(2)) = A(i);
          H(ind(2),ind(1)) = A(i);
       else
          H(ind(1),ind(1)) = A(i); 
       end
    end
end

function b = aux_linPart(A,E)

    % find monomials with polynomial degree 1
    ind = find(sum(E,1) == 1);
    
    % determine slope for each factor
    b = zeros(size(E,1),1);
    
    for i = 1:length(b)
        indTemp = find(E(i,ind) > 0);
        if ~isempty(indTemp)
           b(i) = sum(A(:,ind(indTemp)),2); 
        end
    end
end

function [val,indices] = aux_updateValues(val,indices,val_,ind_)

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
               index = find(all(bsxfun(@withinTol, val_(ind2,i), val(ind1,:)), 1));  
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

function res = aux_subsFactorValue(obj,ind,val)
% substitute factors with a numerical values

    % extract object properties
    E = obj.E;
    G = obj.G;
    EC = obj.EC;
    A = obj.A;

    % substite value into generators
    temp = val.^obj.E(ind,:);
    
    for i = 1:size(temp,1)
        G = G * diag(temp(i,:));
    end
    
    E(ind,:) = [];
    
    % substite value into constraints
    temp = val.^obj.EC(ind,:);
    
    for i = 1:size(temp,1)
        A = A * diag(temp(i,:));
    end
    
    EC(ind,:) = [];
    
    % remove factor identifier
    id = obj.id;
    id(ind) = [];
    
    % construct resulting conPolyZono object
    res = conPolyZono(obj.c,G,E,A,obj.b,EC,obj.GI,id);
    
    % remove redundant monomials
    res = compact_(res,'all',eps);
    
    % remove trivial constriants
    res = aux_removeTrivialConstraints(res);

end

function res = aux_removeTrivialConstraints(obj)
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
        res = conPolyZono(obj.c,obj.G,obj.E,A,b,obj.EC, obj.GI,obj.id);
    else
        res = obj; 
    end
end
    
% ------------------------------ END OF CODE ------------------------------
