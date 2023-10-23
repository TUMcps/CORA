function [Z,varargout] = reduce(Z,method,varargin)
% reduce - reduces the order of a zonotope, the resulting zonotope is an
%    over-approximation of the original zonotope
%
% Syntax:
%    Z = reduce(Z,method,order)
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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%2 inputs
if nargin==2
    order=1;
    filterLength=[];
%3 inputs
elseif nargin==3
    order=varargin{1};
    filterLength=[];
%4 inputs
elseif nargin==4
    order=varargin{1};
    filterLength=varargin{2};
%5 inputs
elseif nargin==5
    order=varargin{1};
    filterLength=varargin{2};
    option = varargin{3};
%6 inputs
elseif nargin==6
    order=varargin{1};
    filterLength=varargin{2};
    option = varargin{3};
    alg = varargin{4};
end

% remove substring necessary for special reduction for polyZonotopes (not
% needed here
if startsWith(method,'approxdep_')
    method = erase(method,'approxdep_');
end

% select option
if strcmp(method,'girard')
    Z=reduceGirard(Z,order);

%option='idx'
elseif strcmp(method,'idx')
    % note: var 'order' is not an order here
    Z = reduceIdx(Z,order);

elseif strcmp(method,'adaptive')
    % note: var 'order' is not an order here!
    [Z,dHerror,gredIdx] = reduceAdaptive(Z,order);
    % additional output arguments
    varargout{1} = dHerror;
    varargout{2} = gredIdx;

elseif strcmp(method,'combastel')
    Z=reduceCombastel(Z,order);

elseif strcmp(method,'pca')
    Z = reducePCA(Z,order);

elseif strcmp(method,'methA')
    Z=reduceMethA(Z,order);

elseif strcmp(method,'methB')
    Z=reduceMethB(Z,order,filterLength); 

elseif strcmp(method,'methC')
    Z=reduceMethC(Z,order,filterLength);

elseif strcmp(method,'methE')
    Z=reduceMethE(Z,order);  

elseif strcmp(method,'methF')
    Z=reduceMethF(Z);   

elseif strcmp(method,'redistribute')
    Z=reduceRedistribute(Z,order);   

elseif strcmp(method,'cluster')
    Z=reduceCluster(Z,order, option);

elseif strcmp(method,'scott')
    Z=reduceScott(Z,order);

elseif strcmp(method,'valero')
    Z=reduceValero(Z,order);

% elseif strcmp(method,'KclusterAllDim')
%     % order must be 1
%     Zred=reduceKclusterAllDim(Z,order);

% elseif strcmp(method,'clusterIter')
%     Zred=reduceClusterIter(Z,order); 

elseif strcmp(method,'constOpt')
    option = 'svd';
    alg = 'interior-point';
    Z = reduceConstOpt(Z,order, option, alg);  

% wrong method
else
    throw(CORAerror('CORA:wrongValue','second',...
        "'adaptive', 'cluster', 'combastel', 'constOpt', 'girard'" + ...
        "'methA', 'methB', 'methC', 'pca', 'scott', 'redistribute', or 'valero'"));
end

% ------------------------------ END OF CODE ------------------------------
