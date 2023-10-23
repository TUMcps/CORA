function Zred = reduceCluster(Z, order, clusterMethod)
% reduceCluster - method to reduce the order of a zonotope
%
% Cluster Methods available:
%    1 - spherical kmeans (matlab implementation)
%    2 - emMean (Cluster after projecting in all d dimensions)
%    3 - emMeanMod (Cluster using kmeans with abs(cos) as sim. measure
%    4 - emLineMean (Clustering using em, line as cluster rep., svd)
%    5 - hierarchicalCluster (Hierarchical clustering using cos. sim., single)
%    6 - hierarchicalCluster (Hierarchical clustering using cos. sim., complete)
%    7 - hierarchicalCluster (Hierarchical clustering using cos. sim., average)
%    8 - hierarchicalCluster (Hierarchical clustering using cos. sim., weighted)
%    9 - kMed (Clustering, using the medoid as center instead of the mean)
%    10 - hybrid methods between Line Clustering (4) and PCA
%    11 - cluster in 2d cluster --> remove half of them
%
% Syntax:
%    Zred = reduceCluster(Z, order, clusterMethod)
%
% Inputs:
%    Z - zonotope object
%    order - order of the reduced zonotope (order=#generators/dim)
%    clusterMethod - see above
%
% Outputs:
%    Zred - reduced zonotope
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Authors:       Anna Kopetzki
% Written:       11-June-2016
% Last update:   27-June-2018
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% initialize Z_red
Zred=Z;

% pick generators to reduce
[center, Gunred, Gred] = pickedGenerators(Z,order);

if ~isempty(Gred)

    % runs of the EmAlgorithm 
    runs=5;
    
    rng(1); % For reproducibility (set seed for kmeans)

    % get dimension dim from zonotope
    d = length(center);
    
    %Delete zero-generators
    G=nonzeroFilter(Gred);

    % robustness: (1,1) and (-1,-1) are the same generator 
    % --> make sure that first dimension !=0 is always >0
    for k=1:size(G,2)
        nonzeroInd=find(G(:,k)); % find all dimension with != 0 
        if G(nonzeroInd(1),k) < 0
            G(:,k)=(-1)*G(:,k);
        end
    end

    % shift the zonotope such that its center is the origin 
    % set cen=origin


    % Cluster the generators in dim Clusters - generator matrix Cgen 
    % Similarity measure: cosine similarity
    % Replicates: how often kmeans is executed
    if clusterMethod == 1
        [ind, gen] = kmeans(G',d, 'Distance', 'cosine');
        Gt = G';
        Cgen=gen';
        for k=1:d
            w = mean(sum(Gt(ind ==k, :).^2, 2));
            Cgen(:,k) = w * Cgen(:,k);
        end
    elseif clusterMethod == 2
        Cgen = emMean(Z, G, d, runs); % cluster all Dim
    elseif clusterMethod == 3
        Cgen = emMeanMod(Z, G, d, runs);
    elseif clusterMethod == 4
        Cgen = emLineMean(Z, G, d, runs);
    elseif clusterMethod == 5
        Cgen = hierarchicalCluster(Z, G, d, runs, clusterMethod);
    elseif clusterMethod == 6
        Cgen = hierarchicalCluster(Z, G, d, runs, clusterMethod);
    elseif clusterMethod == 7
        Cgen = hierarchicalCluster(Z, G, d, runs, clusterMethod);
    elseif clusterMethod == 8
        Cgen = hierarchicalCluster(Z, G, d, runs, clusterMethod);
    elseif clusterMethod == 9
        Cgen = kMed(Z, G, d, runs);
    elseif clusterMethod == 10 % Hybrid of PCA and Line clustering
        Cgen = hybrid(Z, G, d, runs);
    end
    
    % map generators
    Gtrans = Cgen\G;

    % box generators
    Gbox = diag(sum(abs(Gtrans),2));

    % transform generators back
    Gred = Cgen*Gbox;
end

%build reduced zonotope
Zred.c = center;
Zred.G = [Gunred,Gred];


% ------------------------------ END OF CODE ------------------------------
