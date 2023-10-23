function cPZ = reduce(cPZ,method,order)
% reduce - reduces the number of generators of a constrained polynomial 
%    zonotope, resulting in an outer-approximation
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
%    cPZ - reduced conPolyZono object
%
% Example: 
%    c = [0;0];
%    G = [1 0 0.1 0.05 0.01;-2 1 0.2 -0.01 0.01];
%    E = [1 2 3 1 2;0 1 0 2 2;0 0 0 0 0];
%    A = [1 0.25 0.2 0.05 0.01];
%    b = 0.25;
%    EC = [1 0 2 1 2;1 0 0 2 2;0 1 0 0 0];
%    GI = [0 0.05 0.05;0.1 0 0.02];
%    cPZ = conPolyZono(c,G,E,A,b,EC,GI);
%
%    cPZ_ = reduce(cPZ,'girard',6);
%
%    figure; hold on;
%    plot(cPZ,[1,2],'FaceColor','r','Splits',12);
%    plot(cPZ_,[1,2],'b','Splits',12);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/reduce, polyZonotope/reduce

% Authors:       Niklas Kochdumper
% Written:       25-January-2021 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % remove redundant exponents
    cPZ = compact_(cPZ,'all',eps);
    
    % get number of generators, constraints, etc. 
    n = dim(cPZ); h = size(cPZ.G,2); m = size(cPZ.A,1);
    q = size(cPZ.GI,2); z = size(cPZ.A,2);

    % check if reduction is required
    if order*n < h + q + z

        % transform to equivalent higher-dimensional polynomial zonotope
        c = [cPZ.c; -cPZ.b];
        G = blkdiag(cPZ.G,cPZ.A);
        E = [cPZ.E,cPZ.EC];

        GI = cPZ.GI;
        if ~isempty(GI)
            GI = [GI; zeros(m,size(GI,2))];
        end
    
        pZ = polyZonotope(c,G,GI,E,cPZ.id);
        
        % remove redundancies from the representation
        pZ = compact_(pZ,'all',eps);

        % determine required order for high-dimensional poly. zonotope
        order_ = order*n/(2*(n + m));

        % reduce the high-dimensional polynomial zonotope
        pZ = reduce(pZ,method,order_);

        % transform back to constrained polynomial zonotope
        A = []; b = []; E_ = []; GI = [];
        
        if isempty(pZ.G) || all(all(abs(pZ.G(1:n,:)) < eps))
            
            pZ = aux_removeIndepGens(pZ);
            
            c = pZ.c(1:n); G = pZ.G(1:n,:); E = pZ.E; 
            id = (1:length(pZ.id))' + max(pZ.id) + 1;
            if m > 0
               b = -pZ.c(n+1:end);A = pZ.G(n+1:end,:); E_ = pZ.E;
            end
            
        else
            
            c = pZ.c(1:n); G = pZ.G(1:n,:); E = pZ.E;
            GI = pZ.GI; id = pZ.id;

            if m > 0

               % extract constraints
               A = pZ.G(n+1:end,:); b = -pZ.c(n+1:end); E_ = pZ.E;
               Anew = GI(n+1:end,:);

               % try to add uncertainty on constraint coming from the
               % independent generators to the dependent generators
               if strcmp(method,'girard')
                   
                  [A,Anew] = aux_addUncertainty(A,G,pZ.E, ...
                                            Anew(:,end-m+1:end));
                   
                  if ~isempty(Anew)
                      A = [A,Anew];
                      E_ = blkdiag(E_,eye(size(Anew,2)));
                      E = [E;zeros(size(Anew,2),size(E,2))];
                      id = [id; (max(id)+1:max(id)+size(Anew,2))'];              
                  end
                     
               % introduce new dependent factors to represent uncertainty
               % on the constraints    
               else
                   ind = find(sum(abs(Anew),1) > 0);
                   ind_ = setdiff(1:size(Anew,2),ind);
                   Anew = Anew(:,ind);
                   A = [A,Anew];
                   G = [G,GI(1:n,ind)];
                   GI = GI(:,ind_);
                   E_ = blkdiag(E_,eye(size(Anew,2)));
                   E = blkdiag(E,eye(size(Anew,2)));
                   id = [id; (max(id)+1:max(id)+size(Anew,2))'];                    
               end
            end

            if ~isempty(GI)
                GI = GI(1:n,:);
            end
        end
        
        % construct resulting conPolyZono object
        cPZ = conPolyZono(c,G,E,A,b,E_,GI,id);
        
        % remove redundancies
        cPZ = compact_(cPZ,'all',eps);
    end
end


% Auxiliary functions -----------------------------------------------------

function [A,Anew] = aux_addUncertainty(A,G,E,GI)
% try to add uncertainty on constraint coming from the independent 
% generators to the dependent generators

   Anew = []; m = size(A,1); 
   ind = find(sum(E,1) == 1 & sum(abs(G),1) == 0 & sum(abs(A) > 0,1) == 1);

   for i = 1:m
      index = find(abs(A(i,ind)) > 0);
      if ~isempty(index)
          A(i,ind(index(1))) = abs(A(i,ind(index(1)))) + abs(GI(i,i));
      else
          Anew = [Anew, GI(:,i)];
      end
   end
end

function pZ = aux_removeIndepGens(pZ)
% redefine independent generators as new dependent generators

    G = [pZ.G pZ.GI];
    E = blkdiag(pZ.E,eye(size(pZ.GI,2)));
    
    if isempty(pZ.id)
        id = (1:size(pZ.GI,2))';
    else
        temp = max(pZ.id);
        id = [pZ.id; (temp+1:temp+size(pZ.GI,2))'];
    end
    
    pZ = polyZonotope(pZ.c,G,[],E,id);
end

% ------------------------------ END OF CODE ------------------------------
