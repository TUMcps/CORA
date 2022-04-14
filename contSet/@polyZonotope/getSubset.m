function S = getSubset(pZ,id,val)
% getSubset - extract a subset by specifying new ranges for the factors
%
% Syntax:  
%    S = getSubset(pZ,id,val)
%
% Inputs:
%    pZ - polyZonotope object
%    id - vector containing the identifiers of the factors that are changed
%    val - new values for the selected factors (interval or vector) 
%
% Outputs:
%    S - extracted subset (polyZonotope object or single point)
%
% Example: 
%    pZ = polyZonotope([0;0], [2 0 1;1 2 1],[],[1 0 1;0 1 3]);
%    
%    S = getSubset(pZ,2,interval(-0.5,0.9));
%
%    figure
%    hold on
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(S,[1,2],'b','Filled',true,'EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope

% Author:       Niklas Kochdumper
% Written:      23-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % check if value is a vector an interval
    if isnumeric(val)
       
        % loop over all variables that are substituted
        for i = 1:length(id)
           pZ = getSubsetPoint(pZ,id(i),val(i)); 
        end
        
    else
        
        % loop over all variables that are substituted
        for i = 1:length(id)
            if rad(val(i)) ~= 1 || center(val(i)) ~= 0
                pZ = getSubsetInterval(pZ,id(i),val(i)); 
            end
        end
    end
    
    % remove redudancies in the exponent matrix
    [expMat,G] = removeRedundantExponents(pZ.expMat,pZ.G);
    
    % add all constant parts to the center
    ind = find(sum(expMat,1) == 0);
    
    if ~isempty(ind)
        c = pZ.c + sum(G(:,ind),2);
        G(:,ind) = [];
        expMat(:,ind) = [];
    else
        c = pZ.c;
    end
    
    % remove empty rows from the exponent matrix
    id = pZ.id;
    ind = find(sum(expMat,2) == 0);
    if ~isempty(ind)
        expMat(ind,:) = [];
        id(ind) = [];
    end
    
    % construct the resulting subset
    if isempty(G) && isempty(pZ.Grest)
        S = c;
    else
        S = polyZonotope(c,G,pZ.Grest,expMat,id);
    end    
end


% Auxiliary Function ------------------------------------------------------

function pZ = getSubsetPoint(pZ,id,val)

    % find variable that should be substitudet
    ind = find(pZ.id == id);
    
    if isempty(ind)
       error('The specified "id" does not exist!'); 
    end
    
    % modify all generators that are affected by the variable
    for i = 1:size(pZ.expMat,2)
        if pZ.expMat(ind,i) ~= 0
           pZ.G(:,i) = pZ.G(:,i) * val^pZ.expMat(ind,i);
           pZ.expMat(ind,i) = 0;
        end
    end
end

function pZ = getSubsetInterval(pZ,id,val)

    % find variable that should be substitudet
    ind = find(pZ.id == id);
    
    if isempty(ind)
       error('The specified "id" does not exist!'); 
    end
    
    % compute coefficients of all required polynomials (a + b*x)^e
    a = center(val);
    b = rad(val);
    
    expMat = pZ.expMat;
    G = pZ.G;
    
    coeffs = cell(max(expMat(ind,:)),1);
    exps = cell(length(coeffs),1);
    
    coeffs{1} = [b a];
    exps{1} = [1 0];
    
    for i = 2:length(coeffs)
        coeffs{i} = conv(coeffs{i-1},[b,a]);
        exps{i} = i:-1:0;
    end
    
    % modify all generators that are affected by the variable
    for i = 1:size(pZ.expMat,2)
        if pZ.expMat(ind,i) ~= 0
           
           % get coefficients and exponenst of the polnomial
           e = pZ.expMat(ind,i);
           C = coeffs{e};
           E = exps{e};
            
           % consider constant part of the polynomial
           G(:,i) = pZ.G(:,i) * C(end);
           expMat(ind,i) = 0;
           
           % consider remaining part of the polynomial
           expMat_ = pZ.expMat(:,i) * ones(1,length(E)-1);
           expMat_(ind,:) = E(1:end-1);
           
           G_ = pZ.G(:,i) * ones(1,length(E)-1);
           G_ = G_ * diag(C(1:end-1));
           
           expMat = [expMat,expMat_];
           G = [G,G_];
        end
    end
    
    pZ.expMat = expMat;
    pZ.G = G;
end

%------------- END OF CODE --------------