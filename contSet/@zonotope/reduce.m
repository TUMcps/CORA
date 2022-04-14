function Zred = reduce(Z,option,varargin)
% reduce - reduces the order of a zonotope, the resulting zonotope is an
%    over-approximation of the original zonotope
%
% Syntax:  
%    Zred = reduce(Z,option,order)
%
% Inputs:
%    Z - zonotope object
%    option - string specifying the reduction method:
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
%    order - order of reduced zonotope
%
% Outputs:
%    Zred - reduced zonotope (over-approximation)
%
% Example: 
%    Z=zonotope(rand(2,10));
%    figure; hold on;
%    plot(Z,[1,2],'g');
%    Zred=reduce(Z,'girard',2);
%    plot(Zred,[1,2],'r');
%    Zred=reduce(Z,'combastel',2);
%    plot(Zred,[1,2],'b');
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
%
% Other m-files required: none
% Subfunctions: see below
% MAT-files required: none
%
% See also: none

% Author:        Matthias Althoff
% Written:       24-January-2007 
% Last update:   15-September-2007
%                27-June-2018
% Last revision: ---

%------------- BEGIN CODE --------------

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
    method = varargin{3};
%6 inputs
elseif nargin==6
    order=varargin{1};
    filterLength=varargin{2};
    method = varargin{3};
    alg = varargin{4};
end

% remove substring necessary for special reduction for polyZonotopes (not
% needed here
if startsWith(option,'approxdep_')
    option = erase(option,'approxdep_');
end

%option='girard'
if strcmp(option,'girard')
    Zred=reduceGirard(Z,order);
%option='adaptive'
elseif strcmp(option,'adaptive')
    % note: var 'order' is not an order here!
    Zred=reduceAdaptive(Z,order);
%option='combastel'
elseif strcmp(option,'combastel')
    Zred=reduceCombastel(Z,order);
%option='PCA'
elseif strcmp(option,'pca')
    Zred = reducePCA(Z,order);
%option='methA'
elseif strcmp(option,'methA')
    Zred=reduceMethA(Z,order);
%option='methB'
elseif strcmp(option,'methB')
    Zred=reduceMethB(Z,order,filterLength); 
%option='methC'
elseif strcmp(option,'methC')
    Zred=reduceMethC(Z,order,filterLength); 
%option='methE'
elseif strcmp(option,'methE')
    Zred=reduceMethE(Z,order);  
%option='methF'
elseif strcmp(option,'methF')
    Zred=reduceMethF(Z);   
%option='redistribute'
elseif strcmp(option,'redistribute')
    Zred=reduceRedistribute(Z,order);   
% option='cluster'
elseif strcmp(option,'cluster')
    Zred=reduceCluster(Z,order, method);
% option='scott'
elseif strcmp(option,'scott')
    Zred=reduceScott(Z,order);
% % option='KclusterAllDim', order must be 1
% elseif strcmp(option,'KclusterAllDim')
%     Zred=reduceKclusterAllDim(Z,order);
% % option='iter'
% elseif strcmp(option,'clusterIter')
%     Zred=reduceClusterIter(Z,order); 
% option='constOpt'
elseif strcmp(option,'constOpt')
    method = 'svd';
    alg = 'interior-point';
    [Zred]=reduceConstOpt(Z,order, method, alg);  

%wrong argument
else
    error('Invalid reduction method!');
end

%------------- END OF CODE --------------
