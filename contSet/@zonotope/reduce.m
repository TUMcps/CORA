function [Z,varargout] = reduce(Z,method,varargin)
% reduce - reduces the order of a zonotope, the resulting zonotope is an
%    over-approximation of the original zonotope
%
% Syntax:
%    Z = reduce(Z,method)
%    Z = reduce(Z,method,order)
%    Z = reduce(Z,method,order,filterLength)
%    Z = reduce(Z,method,order,filterLength,option)
%    if method from {'minVolume', 'VSF', 'TSF', 'shortestGenVol',
%    'shortestGenTSF', 'scalingDet'}:
%    Z = reduce(Z,method,order,composeVolume)
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
%                   - 'scale'
%                   - 'scaleHausdorff'  
%                   - 'redistribute'
%                   - 'valero'
%                   - 'minVolume'       this and below methods assume
%                   removing one of the original generators and rescaling
%                   the remaining ones leads to the minimal volume (SS)
%                   - 'VSF'
%                   - 'TSF'
%                   - 'shortestGenVol'
%                   - 'shortestGenTSF'
%                   - 'scalingDet'
%    order - order of reduced zonotope
%    filterLength - ???
%    options - ???
%    alg - ???
%    composeVolume - computes volume of over-approximation for minVolume...
%
% Outputs:
%    Z - zonotope object
%    dHerror - (optional, only 'adaptive') over-approximation of the
%              Hausdorff distance between the original and reduced zonotope
%    gredIdx - index of reduced generators
%    minimalVolume - (optional) for methods from {'minVolume', 'VSF', 'TSF',
%                    'shortestGenVol', 'shortestGenTSF', 'scalingDet'}
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
%                11-June-2025 (SS, added further minVolume methods)
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

    case 'scale'
        Z = priv_reduceScale(Z,order);

    case 'scaleHausdorff'
        [Z,dHerror] = priv_reduceScaleHausdorff(Z,order);
        varargout{1} = dHerror;
    
    % case 'KclusterAllDim'
    %     % order must be 1
    %     Zred=reduceKclusterAllDim(Z,order);
    
    % case 'clusterIter'
    %     Zred=reduceClusterIter(Z,order); 
    
    case 'constOpt'
        option = 'svd';
        alg = 'interior-point';
        Z = priv_reduceConstOpt(Z,order, option, alg);  

    case {'minVolume', 'VSF', 'TSF', 'shortestGenVol', 'shortestGenTSF', 'scalingDet'}
        if nargin == 3
            [Z, minimalVolume] = priv_minSubZonotope(Z, 'method', method, 'order', varargin{1});
        else
            [Z, minimalVolume] = priv_minSubZonotope(Z, 'method', method, 'order', varargin{1}, 'composeVolume', varargin{2});
        end
        varargout{1} = minimalVolume;
    
    % wrong method
        otherwise
        throw(CORAerror('CORA:wrongValue','second',...
            "'adaptive', 'adaptive-penven', 'cluster', 'combastel', 'constOpt', 'girard'" + ...
            "'methA', 'methB', 'methC', 'pca', 'scott', 'redistribute', 'sadraddini'" + ...
            "'scale', 'scaleHausdorff', or 'valero'" + ...
            "'minVolume', 'VSF', 'TSF', 'shortestGenVol', 'shortestGenTSF', or 'scalingDet'"));
end

% ------------------------------ END OF CODE ------------------------------
