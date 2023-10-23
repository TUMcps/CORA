function pZ = restructureReduceFull(pZ,order,method)
% restructureReduceFull - Calculates a new representation of a polynomial
%    zonotope through reduction of the independent generators 
%
% Syntax:
%    pZ = restructureReduceFull(pZ,order,method)
%
% Inputs:
%    pZ - polyZonotope object
%    order - desired zonotope order of the dependent factors for the
%            resulting polynomial zonotope 
%    method - reduction technique for linear zonotopes
%             (see zonotope/reduce)
%
% Outputs:
%    pZ - polyZonotope object over-approximating input polynomial zonotope
%
% Example:
%    pZ = polyZonotope([0;0],[1 0 1;1 2 -2],[-1 0.1 -0.5;1.2 0.3 0.2],[1 0 1;0 1 2]);
%    pZnew1 = restructure(pZ,'reduceFullGirard',2);
%    pZnew2 = restructure(pZ,'reduceFullMethC',2);
%
%    figure; hold on;
%    plot(pZnew1,[1,2],'FaceColor','g');
%    plot(pZnew2,[1,2],'FaceColor','b');
%    plot(pZ,[1,2],'FaceColor','r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/restructure, zonotope/reduce

% Authors:       Niklas Kochdumper
% Written:       25-July-2018 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% determine how many generators have to be reduced
dim_x = length(pZ.c);
o = length(pZ.id)/dim_x;
o_ = order - o;

% check if the maximum order is already exceeded
if o_ >= 1
    
    % reduce the independent generators
    zono = zonotope([zeros(dim_x,1),pZ.GI]);
    zono = reduce(zono,method,o_);
    
    % construct the new polynomial zonotope 
    Gzono = generators(zono);
    G = [pZ.G, Gzono];
    E = [pZ.E, zeros(size(pZ.E,1),size(Gzono,2)); ...
              zeros(size(Gzono,2),size(pZ.E,2)), eye(size(Gzono,2))];
    
    pZ = polyZonotope(pZ.c,G,[],E);
    
else
    
    % reduce the independent generators
    zono = zonotope([zeros(dim_x,1),pZ.GI]);
    zono = reduce(zono,method,1);       
    
    % calculate reference zonotope that is added to the generators in
    % order to compare the volumes
    inter = interval(pZ);
    zonoRef = zonotope([zeros(dim_x,1),diag(rad(inter)./100)]);
    
    % calculate the volume for the independent generators
    Gind = generators(zono);
    Vind = zeros(dim_x,1);
    
    for i = 1:size(Gind,2)
        zono_ = zonotope(zonoRef.c,[zonoRef.G,Gind(:,i)]);
        Vind(i) = volume_(interval(zono_)); 
    end
    
    % calculate the volume for the dependent generators
    Vdep = zeros(length(pZ.id),1);
    indicesDep = cell(length(pZ.id),1);
    
    for i = 1:length(Vdep)
        
        % find all generators that that depend on the current factor
        ind = find(pZ.E(i,:) > 0);
        indicesDep{i} = ind';
        pZ_ = polyZonotope(zeros(dim_x,1),pZ.G(:,ind), ...
              generators(zonoRef),pZ.E(:,ind));
        
        zono_ = zonotope(pZ_);
        
        % calculate volume of the zonotope over-approximation
        Vdep(i) = volume_(interval(zono_));
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
    
    pZ_ = polyZonotope(zeros(dim_x,1), pZ.G(:,indDep),[], pZ.E(:,indDep));
    zono_ = zonotope(pZ_);
    
    indTemp = ind2(ind2 > length(pZ.id));
    indTemp = indTemp - length(pZ.id)*ones(size(indTemp));
    
    GI = [zono_.G,Gind(:,indTemp)];
    
    % construct the dependent part of the new polynomial zonotope
    indTemp = setdiff(1:size(pZ.G,2),indDep);
    G1 = pZ.G(:,indTemp);
    E1 = pZ.E(:,indTemp);
    
    indTemp = ind1(ind1 > length(pZ.id));
    indTemp = indTemp - length(pZ.id)*ones(size(indTemp));
    
    G2 = Gind(:,indTemp);
    E2 = eye(size(G2,2));
    
    G = [G1,G2];
    
    if isempty(E1)
        E = E2;
    elseif isempty(E2)
        E = E1;
    else         
        E = blkdiag(E1, E2);
    end 
    
    % remove redundant rows from the exponent matrix
    temp = sum(E,2);
    E = E(temp > 0,:);
    
    % construct the resulting restructured polynomial zonotope
    pZ = polyZonotope(pZ.c+zono_.c, G, GI, E);
end

% ------------------------------ END OF CODE ------------------------------
