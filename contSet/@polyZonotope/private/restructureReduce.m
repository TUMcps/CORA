function pZ = restructureReduce(pZ, order, method, varargin)
% restructureReduce - Calculate a new representation of a polynomial
%    zonotope through reduction of the independent generators 
%
% Syntax:
%    pZ = restructureReduce(pZ, order, method)
%    pZ = restructureReduce(pZ, order, method, genOrder)
%
% Inputs:
%    pZ - polyZonotope object
%    order - desired zonotope order of the dependent factors for the
%            resulting polynomial zonotope 
%    method - reduction technique for linear zonotopes
%             (see zonotope/reduce)
%    genOrder - desired zonotope order of the resulting polynomial zonotope
%
% Outputs:
%    pZ - polyZonotope object over-approximating input polynomial zonotope
%
% Example:
%    pZ = polyZonotope([0;0],[1 0 1;1 2 -2],[-1 0.1 -0.5;1.2 0.3 0.2],[1 0 1;0 1 2]);
%    pZnew1 = restructure(pZ,'reduceGirard',2);
%    pZnew2 = restructure(pZ,'reduceMethC',2);
%
%    hold on
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

% parse input arguments
genOrder = setDefaultValues({Inf},varargin);

% check if the maximum zonotope order is exceeded
dim_x = length(pZ.c);
o = (length(pZ.id)+dim_x)/dim_x;

if o <= order    % max order satisfied
    
    % check if additional generators need to be removed
    o_ = (size(pZ.G,2) + dim_x) - genOrder*dim_x;
    
    if o_ > 0
        
        % half the generator length for exponents that are all even
        Gtemp = pZ.G;         
        temp = prod(ones(size(pZ.E))-mod(pZ.E,2),1);
        ind = find(temp == 1);
        Gtemp(:,ind) = 0.5 * Gtemp(:,ind);
       
        % determine length of the generators
        len = sum(Gtemp.^2,1);
        [~,ind] = sort(len,'ascend');
        
        % reduce the smallest generators
        ind = ind(1:o_);
        
        Grem = pZ.G(:,ind);
        ERem = pZ.E(:,ind);
        
        pZ.G(:,ind) = [];
        pZ.E(:,ind) = [];  
        
        % reduce the polynomial zonotope that corresponds to the
        % generators that are removed
        pZ_ = polyZonotope(zeros(dim_x,1),Grem,pZ.GI,ERem);
        
        zono = zonotope(pZ_);
        zono = reduce(zono,method,1);
        
    else       
    
        % reduce the zonotope that corresponds to the independent generators
        zono_ = zonotope([zeros(dim_x,1),pZ.GI]);
        zono = reduce(zono_,method,1);
    end
    
    
    % construct the restructured polynomial zonotope
    Gzono = zono.G;
    G = [pZ.G, Gzono];
    E = [pZ.E, zeros(size(pZ.E,1),size(Gzono,2)); ...
              zeros(size(Gzono,2),size(pZ.E,2)), eye(size(Gzono,2))];

    pZ = polyZonotope(pZ.c+center(zono),G,[],E);
    
else            % max order exceeded
    
    % number of dependent generators that need to be removed
    n = ceil(length(pZ.id) + dim_x - dim_x*order); 
    
    % calculate reference zonotope that is added to the generators in
    % order to compare the volumes
    inter = interval(pZ);
    zonoRef = zonotope([zeros(dim_x,1),diag(rad(inter)./100)]);
    
    % calculate the volume for all dependent generators
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
    [~,ind] = sort(Vdep,'ascend');
    
    ind1 = ind(1:n);
    
    % determine the indices of all generators that are removed
    indicesDep_ = indicesDep(ind1);
    indDep = unique(vertcat(indicesDep_{:}));
    
    Grem = pZ.G(:,indDep);
    pZ.G(:,indDep) = [];
    
    ERem = pZ.E(:,indDep);
    pZ.E(:,indDep) = [];
    
    % check if additional generators need to be removed
    o_ = (size(pZ.G,2) + dim_x) - genOrder*dim_x;
    
    if o_ > 0
       
        % half the generator length for exponents that are all even
        Gtemp = pZ.G;         
        temp = prod(ones(size(pZ.E))-mod(pZ.E,2),1);
        ind = find(temp == 1);
        Gtemp(:,ind) = 0.5 * Gtemp(:,ind);
        
        % determine length of the generators
        len = sum(Gtemp.^2,1);
        [~,ind] = sort(len,'ascend');
        
        % reduce the smallest generators
        ind = ind(1:o_);
        
        Grem = [Grem, pZ.G(:,ind)];
        ERem = [ERem, pZ.E(:,ind)];
        
        pZ.G(:,ind) = [];
        pZ.E(:,ind) = [];   
    end
    
    % construct a polynomial zonotope with all generators that are
    % removed and reduce it to order 1
    pZ_ = polyZonotope(zeros(dim_x,1), Grem ,pZ.GI, ERem);
    zono_ = zonotope(pZ_);
    zono = reduce(zono_,method,1);
    
    if ~isempty(pZ.E)
        pZ.E(ind1,:) = [];
    end
    
    % construct the restructured polynomial zonotope
    Gzono = zono.G;
    G = [pZ.G, Gzono];
    E = blkdiag(pZ.E, eye(size(Gzono,2)));
    
    pZ = polyZonotope(pZ.c+zono.c,G,[],E);   
end

% ------------------------------ END OF CODE ------------------------------
