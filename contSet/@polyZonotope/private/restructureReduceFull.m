function res = restructureReduceFull(pZ, order, method)
% restructureReduce - Calculate a new representation of a polynomial
%                     zonotope through reduction of the independent 
%                     generators 
%
% Syntax:  
%    res = restructureReduce(pZ, order, method)
%
% Inputs:
%    pZ - polyZonotope object
%    order - desired zonotope order of the dependent factors for the
%            resulting polynomial zonotope 
%    method - reduction technique for linear zonotopes
%             (see zonotope/reduce)
%
% Outputs:
%    res - resuling polyZonotope object which over-approximates pZ
%
% Example:
%    pZ = polyZonotope([0;0],[1 0 1;1 2 -2],[-1 0.1 -0.5;1.2 0.3 0.2],[1 0 1;0 1 2]);
%    pZnew1 = restructure(pZ,'reduceFullGirard',2);
%    pZnew2 = restructure(pZ,'reduceFullMethC',2);
%
%    hold on
%    plot(pZnew1,[1,2],'g','Filled',true,'EdgeColor','none');
%    plot(pZnew2,[1,2],'b','Filled',true,'EdgeColor','none');
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/restructure, zonotope/reduce

% Author:       Niklas Kochdumper
% Written:      25-July-2018 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % determine how many generators have to be reduced
    dim_x = length(pZ.c);
    o = length(pZ.id)/dim_x;
    o_ = order - o;

    % check if the maximum order is already exceeded
    if o_ >= 1
        
        % reduce the independent generators
        zono = zonotope([zeros(dim_x,1),pZ.Grest]);
        zono = reduce(zono,method,o_);
        
        % construct the new polynomial zonotope 
        Gzono = generators(zono);
        G = [pZ.G, Gzono];
        expMat = [pZ.expMat, zeros(size(pZ.expMat,1),size(Gzono,2)); ...
                  zeros(size(Gzono,2),size(pZ.expMat,2)), eye(size(Gzono,2))];
        
        res = polyZonotope(pZ.c,G,[],expMat);
        
    else
        
        % reduce the independent generators
        zono = zonotope([zeros(dim_x,1),pZ.Grest]);
        zono = reduce(zono,method,1);       
        
        % calculate reference zonotope that is added to the generators in
        % order to compare the volumes
        inter = interval(pZ);
        zonoRef = zonotope([zeros(dim_x,1),diag(rad(inter)./100)]);
        
        % calculate the volume for the independent generators
        Gind = generators(zono);
        Vind = zeros(dim_x,1);
        
        for i = 1:size(Gind,2)
            zono_ = zonotope([zonoRef.Z,Gind(:,i)]);
            Vind(i) = volume(interval(zono_)); 
        end
        
        % calculate the volume for the dependent generators
        Vdep = zeros(length(pZ.id),1);
        indicesDep = cell(length(pZ.id),1);
        
        for i = 1:length(Vdep)
            
            % find all generators that that depend on the current factor
            ind = find(pZ.expMat(i,:) > 0);
            indicesDep{i} = ind';
            pZ_ = polyZonotope(zeros(dim_x,1),pZ.G(:,ind), ...
                  generators(zonoRef),pZ.expMat(:,ind));
            
            zono_ = zonotope(pZ_);
            
            % calculate volume of the zonotope over-approximation
            Vdep(i) = volume(interval(zono_));
        end
        
        % find generators with the smallest volume => smallest
        % over-approximation by removal
        V = [Vdep;Vind];
        [~,ind] = sort(V,'descend');
        
        ind1 = ind(1:order*dim_x);
        ind2 = ind(order*dim_x + 1:end);
        
        
        % construct the independent part of the new polynomial zonotope
        indTemp = ind2(ind2 <= length(pZ.id));
        indicesDep_ = indicesDep(indTemp);
        indDep = unique(vertcat(indicesDep_{:}));
        
        pZ_ = polyZonotope(zeros(dim_x,1), pZ.G(:,indDep),[], pZ.expMat(:,indDep));
        zono_ = zonotope(pZ_);
        
        indTemp = ind2(ind2 > length(pZ.id));
        indTemp = indTemp - length(pZ.id)*ones(size(indTemp));
        
        c = center(zono_);
        Grest = [generators(zono_),Gind(:,indTemp)];
        
        % construct the dependent part of the new polynomial zonotope
        indTemp = setdiff(1:size(pZ.G,2),indDep);
        G1 = pZ.G(:,indTemp);
        expMat1 = pZ.expMat(:,indTemp);
        
        indTemp = ind1(ind1 > length(pZ.id));
        indTemp = indTemp - length(pZ.id)*ones(size(indTemp));
        
        G2 = Gind(:,indTemp);
        expMat2 = eye(size(G2,2));
        
        G = [G1,G2];
        
        if isempty(expMat1)
            expMat = expMat2;
        elseif isempty(expMat2)
            expMat = expMat1;
        else         
            expMat = [expMat1, zeros(size(expMat1,1),size(expMat2,2)); ...
                      zeros(size(expMat2,1),size(expMat1,2)), expMat2 ];
        end 
        
        % remove redundant rows from the exponent matrix
        temp = sum(expMat,2);
        expMat = expMat(temp > 0,:);
        
        % construct the resulting restructured polynomial zonotope
        res = polyZonotope(pZ.c+c, G, Grest, expMat);
    end

%------------- END OF CODE --------------