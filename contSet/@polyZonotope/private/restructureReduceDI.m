function pZ = restructureReduceDI(pZ, dfOrder, method, varargin)
% restructureReduceDI - Calculate a new representation of a polynomial
%    zonotope through simultanious reduction of dependent and independent
%    generators
%
% Syntax:
%    pZ = restructureReduceDI(pZ, order, method)
%    pZ = restructureReduceDI(pZ, order, method, genOrder)
%
% Inputs:
%    pZ - polyZonotope object
%    dfOrder - desired zonotope order of the dependent factors for the
%            resulting polynomial zonotope 
%    method - reduction technique for linear zonotopes
%             (see zonotope/reduce)
%    genOrder - (optional) desired zonotope order of the resulting
%               polynomial zonotope
%
% Outputs:
%    pZ - polyZonotope object over-approximating input polynomial zonotope
%
% Example:
%    pZ = polyZonotope([0;0],[1 0 1;1 2 -2],[-1 0.1 -0.5;1.2 0.3 0.2],[1 0 1;0 1 2]);
%    pZnew1 = restructure(pZ,'reduceDIGirard',2);
%    pZnew2 = restructure(pZ,'reduceDIMethC',2);
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

% Authors:       Victor Gassmann
% Written:       29-October-2021 
% Last update:   26-October-2022 (VG, retain existing ids)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    genOrder = setDefaultValues({Inf},varargin);
    
    % make sure that polyZonotope is irreducible
    pZ = compact_(pZ,'states',eps);
    
    % read out dimension of polynomial zonotope
    n = dim(pZ);
    
    % construct equivalent polyZonotope without indep. generators
    G = [pZ.G,pZ.GI];
    E = blkdiag(pZ.E,eye(size(pZ.GI,2)));
    if ~isempty(pZ.id)
        id = [pZ.id;max(abs(pZ.id))+(1:size(pZ.GI,2))'];
    else
        id = (1:size(pZ.GI,2))';
    end
    [d,m] = size(E);
    
    % case 1: everything is within limits
    if d <= dfOrder*n && m <= genOrder*n
        pZ.G = G;
        pZ.GI = zeros(n,0);
        pZ.E = E;
        pZ.id = id;
        return;
    end
    
    % either too many dep. factors, or too many generators, or both
    
    
    % case 2: too many generators, but still fit inside dfOrder budget
    if m>genOrder*n && d<=(dfOrder-1)*n
        
        % number of generators that need to be reduced to n generators
        M = m-(genOrder-1)*n;
        
        ind_r = aux_gens2reduce(G,E,M);
        
        % reduce resulting generators to n generators
        G_ = generators(reduce(zonotope(zeros(n,1),G(:,ind_r)),method,1));
        m = size(G_,2);
        pZ.G = [G(:,~ind_r),G_];
        pZ.GI = zeros(n,0);
        pZ.E = blkdiag(E(:,~ind_r),eye(m));
        pZ.id = [id;max(abs(id))+(1:m)'];
        pZ = compact_(pZ,'states',eps);
                        
    %case 3+4: generators fit/dont fit, but too many dep. factors
    elseif d>dfOrder*n
        
        % case 3+4
        D = d-(dfOrder-1)*n;
        
        % calculate the volume for all dependent generators
        V = zeros(d,1);
        Ind = false(d,m);
        for i = 1:d
            
            % find all generators that depend on the current factor
            ind = E(i,:) > 0;
            Ind(i,:) = ind;
            pZ_ = polyZonotope(zeros(n,1),G(:,ind), ...
                  zeros(n,0),E(:,ind),id);
            
            % calculate "volume" of the zonotope over-approximation
            V(i) = sum(rad(interval(zonotope(pZ_))));
        end
        
        % find generators with the smallest volume => smallest
        % over-approximation by removal
        [~,ii] = sort(V,'ascend');
        
        % select all generators that need to be removed
        indMask = any(Ind(ii(1:D),:),1);
        
        % number of generators currently (if we were to construct res now)
        m_ = sum(~indMask)+n;
        
        % case 4: too many generators
        if m_>genOrder*n
            % remove further generators
            % generators to remove additionally
            N = m_-(genOrder-1)*n;
            
            tmp = 1:m;
            ii_mask_not = tmp(~indMask);
            
            % find indices of generators to reduce
            ind_r = aux_gens2reduce(G(:,~indMask),E(:,~indMask),N);
            
            % convert ind_r to mask in original "G" space
            ind_gr = ismember(1:m,ii_mask_not(ind_r));
            
            % mask of all generators to reduce
            indMask = indMask & ind_gr;
        end    
        
        Zr = reduce(zonotope(polyZonotope(zeros(n,1),G(:,indMask),...
                    zeros(n,0),E(:,indMask))),method,1);
        Gr = generators(Zr);
        m = size(Gr,2);
        G = [G(:,~indMask),Gr];
        E = blkdiag(E(ii(D+1:end),~indMask),eye(m));
        id = [id(ii(D+1:end));
              max(abs(id(ii(D+1:end))))+(1:m)'];
        % construct resulting polyZonotope
        pZ = polyZonotope(pZ.c+Zr.c, G, zeros(n,0), E, id);
        pZ = compact_(pZ,'states',eps);
    else
        throw(CORAerror('CORA:specialError',...
            'Bug in restructureReduceDI: Check logic!'));
    end
end


% Auxiliary functions -----------------------------------------------------

function ind_r = aux_gens2reduce(G,E,M)
    % sort generators according to metric
    m = size(G,2);
    Gtemp = G;         
    ind = all(mod(E,2)==0,1);
    Gtemp(:,ind) = 0.5 * Gtemp(:,ind);
    
    % determine length of the generators
    len = sum(Gtemp.^2,1);
    [~,ii_s] = sort(len,'ascend');
    ind_r = ismember(1:m,ii_s(1:M));
end

% ------------------------------ END OF CODE ------------------------------
