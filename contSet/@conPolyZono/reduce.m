function cPZ = reduce(cPZ,method,order)
% reduce - reduces the number of generators of a constrained polynomial 
%          zonotope
%
% Syntax:  
%    cPZ = reduce(cPZ,method,order)
%
% Inputs:
%    cPZ - conPolyZono object
%    method - reduction algorithm (see zonotope/reduce)
%    order - desired order of reduced constrained polynomial zonotope
%
% Outputs:
%    pZ - reduced polynomial zonotope
%
% Example: 
%    c = [0;0];
%    G = [1 0 0.1 0.05 0.01;-2 1 0.2 -0.01 0.01];
%    expMat = [1 2 3 1 2;0 1 0 2 2;0 0 0 0 0];
%    A = [1 0.25 0.2 0.05 0.01];
%    b = 0.25;
%    expMat_ = [1 0 2 1 2;1 0 0 2 2;0 1 0 0 0];
%    Grest = [0 0.05 0.05;0.1 0 0.02];
%    cPZ = conPolyZono(c,G,expMat,A,b,expMat_,Grest);
%
%    cPZ_ = reduce(cPZ,'girard',6);
%
%    figure; hold on;
%    plot(cPZ,[1,2],'r','Filled',true,'Splits',15);
%    plot(cPZ_,[1,2],'b','Splits',15);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/reduce, polyZonotope/reduce

% Author:       Niklas Kochdumper
% Written:      25-January-2021 
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % remove redundant exponents
    cPZ = compact(cPZ);
    
    % get number of generators, constraints, etc. 
    n = dim(cPZ); h = size(cPZ.G,2); m = size(cPZ.A,1);
    q = size(cPZ.Grest,2); z = size(cPZ.A,2);

    % check if reduction is required
    if order*n < h + q + z

        % transform to equivalent higher-dimensional polynomial zonotope
        c = [cPZ.c; -cPZ.b];
        G = blkdiag(cPZ.G,cPZ.A);
        expMat = [cPZ.expMat,cPZ.expMat_];

        Grest = cPZ.Grest;
        if ~isempty(Grest)
            Grest = [Grest; zeros(m,size(Grest,2))];
        end
    
        pZ = polyZonotope(c,G,Grest,expMat,cPZ.id);
        
        % remove redundancies from the representation
        pZ = compact(pZ);

        % determine required order for high-dimensional poly. zonotope
        order_ = order*n/(2*(n + m));

        % reduce the high-dimensional polynomial zonotope
        pZ = reduce(pZ,method,order_);

        % transform back to constrained polynomial zonotope
        A = []; b = []; expMat_ = []; Grest = [];
        
        if isempty(pZ.G) || all(all(abs(pZ.G(1:n,:) < eps)))
            
            pZ = removeIndepGens(pZ);
            
            c = pZ.c(1:n); G = pZ.G(1:n,:); expMat = pZ.expMat; 
            id = (1:length(pZ.id))' + max(pZ.id) + 1;
            if m > 0
               b = -pZ.c(n+1:end);A = pZ.G(n+1:end,:); expMat_ = pZ.expMat;
            end
            
        else
            
            c = pZ.c(1:n); G = pZ.G(1:n,:); expMat = pZ.expMat;
            Grest = pZ.Grest; id = pZ.id;

            if m > 0

               % extract constraints
               A = pZ.G(n+1:end,:); b = -pZ.c(n+1:end); expMat_ = pZ.expMat;
               Anew = Grest(n+1:end,:);

               % try to add uncertainty on constraint coming from the
               % independent generators to the dependent generators
               if strcmp(method,'girard')
                   
                  [A,Anew] = addUncertainty(A,G,pZ.expMat, ...
                                            Anew(:,end-m+1:end));
                   
                  if ~isempty(Anew)
                      A = [A,Anew];
                      expMat_ = blkdiag(expMat_,eye(size(Anew,2)));
                      expMat = [expMat;zeros(size(Anew,2),size(expMat,2))];
                      id = [id; (max(id)+1:max(id)+size(Anew,2))'];              
                  end
                     
               % introduce new dependent factors to represent uncertainty
               % on the constraints    
               else
                   ind = find(sum(abs(Anew),1) > 0);
                   ind_ = setdiff(1:size(Anew,2),ind);
                   Anew = Anew(:,ind);
                   A = [A,Anew];
                   G = [G,Grest(1:n,ind)];
                   Grest = Grest(:,ind_);
                   expMat_ = blkdiag(expMat_,eye(size(Anew,2)));
                   expMat = blkdiag(expMat,eye(size(Anew,2)));
                   id = [id; (max(id)+1:max(id)+size(Anew,2))'];                    
               end
            end

            if ~isempty(Grest)
                Grest = Grest(1:n,:);
            end
        end
        
        % construct resulting conPolyZono object
        cPZ = conPolyZono(c,G,expMat,A,b,expMat_,Grest,id);
        
        % remove redundancies
        cPZ = compact(cPZ);
    end
end


% Auxiliary Functions -----------------------------------------------------

function [A,Anew] = addUncertainty(A,G,expMat,Grest)
% try to add uncertainty on constraint coming from the independent 
% generators to the dependent generators

   Anew = []; m = size(A,1); 
   ind = find(sum(expMat,1) == 1 & sum(abs(G),1) == 0 & sum(abs(A) > 0,1) == 1);

   for i = 1:m
      index = find(abs(A(i,ind)) > 0);
      if ~isempty(index)
          A(i,ind(index(1))) = abs(A(i,ind(index(1)))) + abs(Grest(i,i));
      else
          Anew = [Anew, Grest(:,i)];
      end
   end
end

function pZ = removeIndepGens(pZ)
% redefine independent generators as new dependent generators

    G = [pZ.G pZ.Grest];
    expMat = blkdiag(pZ.expMat,eye(size(pZ.Grest,2)));
    
    if isempty(pZ.id)
        id = (1:size(pZ.Grest,2))';
    else
        temp = max(pZ.id);
        id = [pZ.id; (temp+1:temp+size(pZ.Grest,2))'];
    end
    
    pZ = polyZonotope(pZ.c,G,[],expMat,id);
end

%------------- END OF CODE --------------  