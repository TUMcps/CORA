function [Z,varargout] = reduce(Z,method,varargin)
% reduce - reduces the order of a zonotope, the resulting zonotope is an
%    over-approximation of the original zonotope
%
% Syntax:
%    Z = reduce(Z,method)
%    Z = reduce(Z,method,order)
%    Z = reduce(Z,method,order,filterLength)
%    Z = reduce(Z,method,order,filterLength,option)
%
% Inputs:
%    Z - zonotope object
%    method - string specifying the reduction method:
%                   - 'adaptive'        Thm. 3.2. in [6]
%                   - 'cluster'         Sec. III.B in [3]
%                   - 'combastel'       Sec. 3.2 in [4]
%                   - 'constOpt'        Sec. III.D in [3]
%                   - 'girard'          Sec. 4 in [2]
%                   - 'methA'           Sec. 2.5.5 in [1]
%                   - 'methB'           Sec. 2.5.5 in [1]
%                   - 'methC'           Sec. 2.5.5 in [1]
%                   - 'pca'             Sec. III.A in [3]
%                   - 'scott'           Appendix of [5]
%                   - 'sadraddini'      Proposition 6 in [8]
%                   - 'redistribute'
%                   - 'valero'          
%    order - order of reduced zonotope
%    filterLength - ???
%    options - ???
%    alg - ???
%
% Outputs:
%    Z - zonotope object
%    dHerror - (optional, only 'adaptive') over-approximation of the
%              Hausdorff distance between the original and reduced zonotope
%    gredIdx - index of reduced generators
%
% Example: 
%    Z = zonotope([1;-1],[2 3 1 -1 -2 1 -4 3 0; 2 3 -2 1 -3 2 1 0 2]);
%    Zred1 = reduce(Z,'girard',2);
%    Zred2 = reduce(Z,'combastel',2);
%
%    figure; hold on;
%    plot(Z,[1,2],'g');
%    plot(Zred1,[1,2],'r');
%    plot(Zred2,[1,2],'b');
%
% References:
%    [1] M. Althoff. "Reachability analysis and its application to the 
%        safety assessment of autonomous cars", 2010
%    [2] A. Girard. "Reachability of uncertain linear systems using
%        zonotopes". 2005
%    [3] A. Kopetzki et al. "Methods for order reduction of zonotopes", 
%        CDC 2017
%    [4] C. Combastel. "A state bounding observer based on zonotopes",
%        ECC 2003
%    [5] J. Scott et al. "Constrained zonotopes: A new tool for set-based 
%        estimation and fault detection", Automatica 2016
%    [6] M. Wetzlinger et al. "Adaptive parameter tuning for reachability
%        analysis of nonlinear systems", HSCC 2021.
%    [7] C.E. Valero et al. "On minimal volume zonotope order reduction",
%        Automatica 2021 (in revision)
%    [8] S. Sadraddini and R. Tedrake. "Linear Encodings for Polytope
%        Containment Problems", CDC 2019 (ArXiV version)
%
% Other m-files required: none
% Subfunctions: see below
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       24-January-2007 
% Last update:   15-September-2007
%                27-June-2018
% Last revision: 06-October-2024 (MW, refactor including priv_)

% ------------------------------ BEGIN CODE -------------------------------

narginchk(2,5);
[order,filterLength,option] = setDefaultValues({1,[],[]},varargin);

% remove substring necessary for special reduction for polyZonotopes (not
% needed here
if startsWith(method,'approxdep_')
    method = erase(method,'approxdep_');
end
% select option
switch method
    case 'girard'
        Z = priv_reduceGirard(Z,order);

    case 'idx'
        % note: var 'order' is not an order here
        Z = priv_reduceIdx(Z,order);

    case 'adaptive'
        % note: var 'order' is not an order here!
        [Z,dHerror,gredIdx] = priv_reduceAdaptive(Z,order);
        % additional output arguments
        varargout{1} = dHerror;
        varargout{2} = gredIdx;

    case strcmp(method,'adaptive-penven')
        % note: var 'order' is not an order here!
        [Z,dHerror,gredIdx] = reduceAdaptive(Z,order,'penven');
        % additional output arguments
        varargout{1} = dHerror;
        varargout{2} = gredIdx;
    
    case 'combastel'
        Z = priv_reduceCombastel(Z,order);
    
    case 'pca'
        Z = priv_reducePCA(Z,order);
    
        % methX
    case 'methA'
        Z = priv_reduceMethA(Z,order);
    
    case 'methB'
        Z = priv_reduceMethB(Z,order,filterLength); 
    
    case 'methC'
        Z = priv_reduceMethC(Z,order,filterLength);
    
    case 'methE'
        Z = priv_reduceMethE(Z,order);  
    
    case 'methF'
        Z = priv_reduceMethF(Z);   
    
        % ---
    case 'redistribute'
        Z = priv_reduceRedistribute(Z,order);   
    
    case 'cluster'
        Z = priv_reduceCluster(Z,order, option);
    
    case 'scott'
        Z = priv_reduceScott(Z,order);
    
    case 'valero'
        Z = priv_reduceValero(Z,order);

    case 'sadraddini'
        Z = priv_reduceSadraddini(Z,order);
    
    % case 'KclusterAllDim'
    %     % order must be 1
    %     Zred=reduceKclusterAllDim(Z,order);
    
    % case 'clusterIter'
    %     Zred=reduceClusterIter(Z,order); 
    
    case 'constOpt'
        option = 'svd';
        alg = 'interior-point';
        Z = priv_reduceConstOpt(Z,order, option, alg);  
    
    % wrong method
        otherwise
        throw(CORAerror('CORA:wrongValue','second',...
            "'adaptive', 'adaptive-penven', 'cluster', 'combastel', 'constOpt', 'girard'" + ...
            "'methA', 'methB', 'methC', 'pca', 'scott', 'redistribute', 'sadraddini', or 'valero'"));
end

% ------------------------------ END OF CODE ------------------------------
