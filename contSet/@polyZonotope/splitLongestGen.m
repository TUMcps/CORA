function [pZsplit,factor] = splitLongestGen(pZ,varargin)
% splitLongestGen - Splits the longest generator dependent generator with a 
%                   polynomial order of 1 for a polynomial zonotope
%
% Syntax:  
%    [pZsplit,factor] = splitLongestGen(pZ)
%    [pZsplit,factor] = splitLongestGen(pZ,polyOrd)
%
% Inputs:
%    pZ - polyZonotope object
%    polyOrd - maximum number of polynomial terms that are splitted exactly
%              (without an over-approximation)
%
% Outputs:
%    pZsplit - cell array of split polyZonotopes
%    factro - identifier of the dependent factor that is split
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%
%    temp = splitLongestGen(pZ);
%    pZsplit1 = splitLongestGen(temp{1});
%    pZsplit2 = splitLongestGen(temp{2});
%
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%
%    figure
%    hold on
%    plot(pZsplit1{1},[1,2],'b','Filled',true,'EdgeColor','none');
%    plot(pZsplit1{2},[1,2],'g','Filled',true,'EdgeColor','none');
%    plot(pZsplit2{1},[1,2],'m','Filled',true,'EdgeColor','none');
%    plot(pZsplit2{2},[1,2],'c','Filled',true,'EdgeColor','none');
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: split, splitDepFactor

% Author:       Niklas Kochdumper
% Written:      29-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % determine longest generator
    len = sum(pZ.G.^2,1);
    [~,ind] = max(len);

    % find factor with the largest exponent
    [~,index] = max(pZ.expMat(:,ind));
    index = pZ.id(index);

    % split the zonotope at the determined generator
    if nargin == 2
       pZsplit = splitDepFactor(pZ,index,varargin{1}); 
    else
       pZsplit = splitDepFactor(pZ,index); 
    end
    
    factor = index;

%------------- END OF CODE --------------