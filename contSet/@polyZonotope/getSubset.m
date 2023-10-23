function S = getSubset(pZ,id,val)
% getSubset - extracts a subset by specifying new ranges for the factors
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
%    figure; hold on;
%    plot(pZ,[1,2],'FaceColor','r');
%    plot(S,[1,2],'FaceColor','b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope

% Authors:       Niklas Kochdumper
% Written:       23-March-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % check input arguments
    inputArgsCheck({{pZ,'att','polyZonotope'};
                    {id,'att','numeric','vector'};
                    {val,'att',{'numeric','interval'}}});
    if representsa_(val,'emptySet',eps)
        throw(CORAerror("CORA:wrongValue","third",'Should not be empty.'))
    end

    % check if value is a vector an interval
    if isnumeric(val)
       
        % loop over all variables that are substituted
        for i = 1:length(id)
           pZ = aux_getSubsetPoint(pZ,id(i),val(i)); 
        end
        
    else
        
        % loop over all variables that are substituted
        for i = 1:length(id)
            if rad(val(i)) ~= 1 || center(val(i)) ~= 0
                pZ = aux_getSubsetInterval(pZ,id(i),val(i)); 
            end
        end
    end
    
    % remove redudancies in the exponent matrix
    [E,G] = removeRedundantExponents(pZ.E,pZ.G);
    
    % add all constant parts to the center
    ind = find(sum(E,1) == 0);
    
    if ~isempty(ind)
        c = pZ.c + sum(G(:,ind),2);
        G(:,ind) = [];
        E(:,ind) = [];
    else
        c = pZ.c;
    end
    
    % remove empty rows from the exponent matrix
    id = pZ.id;
    ind = find(sum(E,2) == 0);
    if ~isempty(ind)
        E(ind,:) = [];
        id(ind) = [];
    end
    
    % construct the resulting subset
    if isempty(G) && isempty(pZ.GI)
        S = c;
    else
        S = polyZonotope(c,G,pZ.GI,E,id);
    end    
end


% Auxiliary functions -----------------------------------------------------

function pZ = aux_getSubsetPoint(pZ,id,val)

    % find variable that should be substitudet
    ind = find(pZ.id == id);
    
    if isempty(ind)
        throw(CORAerror('CORA:wrongValue','second',...
            'be a vector containing the identifiers of the factors that are changed'));
    end
    
    % modify all generators that are affected by the variable
    for i = 1:size(pZ.E,2)
        if pZ.E(ind,i) ~= 0
           pZ.G(:,i) = pZ.G(:,i) * val^pZ.E(ind,i);
           pZ.E(ind,i) = 0;
        end
    end

end

function pZ = aux_getSubsetInterval(pZ,id,val)

    % find variable that should be substitudet
    ind = find(pZ.id == id);
    
    if isempty(ind)
        throw(CORAerror('CORA:wrongValue','second',...
            'be a vector containing the identifiers of the factors that are changed'));
    end
    
    % compute coefficients of all required polynomials (a + b*x)^e
    a = center(val);
    b = rad(val);
    
    E = pZ.E;
    G = pZ.G;
    
    coeffs = cell(max(E(ind,:)),1);
    exps = cell(length(coeffs),1);
    
    coeffs{1} = [b a];
    exps{1} = [1 0];
    
    for i = 2:length(coeffs)
        coeffs{i} = conv(coeffs{i-1},[b,a]);
        exps{i} = i:-1:0;
    end
    
    % modify all generators that are affected by the variable
    for i = 1:size(pZ.E,2)
        if pZ.E(ind,i) ~= 0
           
           % get coefficients and exponenst of the polnomial
           e = pZ.E(ind,i);
           C = coeffs{e};
           Ei = exps{e};
            
           % consider constant part of the polynomial
           G(:,i) = pZ.G(:,i) * C(end);
           E(ind,i) = 0;
           
           % consider remaining part of the polynomial
           E_ = pZ.E(:,i) * ones(1,length(Ei)-1);
           E_(ind,:) = Ei(1:end-1);
           
           G_ = pZ.G(:,i) * ones(1,length(Ei)-1);
           G_ = G_ * diag(C(1:end-1));
           
           E = [E,E_];
           G = [G,G_];
        end
    end
    
    pZ.E = E;
    pZ.G = G;

end

% ------------------------------ END OF CODE ------------------------------
