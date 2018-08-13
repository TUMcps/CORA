function [Zred]=reduce(Z,option,varargin)
% reduce - Reduces the order of a zonotope
% options: 
% cluster 
% combastel
% constOpt
% girard
% methA
% methB
% methC
% pca
% redistribute
% scott
%
% Syntax:  
%    [Zred]=reduce(Z,option,order)
%
% Inputs:
%    Z - zonotope object
%    option - 'girard' or 'althoff'
%    order - order of reduced zonotope
%
% Outputs:
%    Zred - reduced zonotope
%
% Example: 
%    Z=zonotope(rand(2,10));
%    plot(Z,[1,2],'g');
%    hold on
%    Zred=reduce(Z,'girard',2);
%    plot(Zred,[1,2],'r');
%    Zred=reduce(Z,'combastel',2);
%    plot(Zred,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: see below
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      24-January-2007 
% Last update:  15-September-2007
%               27-June-2018
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


%option='girard'
if strcmp(option,'girard')
    Zred=reduceGirard(Z,order);
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
    disp('Error: Second argument is wrong');
    Zred=[];
end

%------------- END OF CODE --------------
