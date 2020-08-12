function [pZsplit] = splitLongestGen(pZ,varargin)
% splitLongestGen - Splits the longest generator dependent generator with a 
%                   polynomial order of 1 for a polynomial zonotope
%
% Syntax:  
%    [pZsplit] = splitOneGen(pZ)
%    [pZsplit] = splitOneGen(pZ,polyOrd)
%
% Inputs:
%    pZ - polyZonotope object
%    polyOrd - maximum number of polynomial terms that are splitted exactly
%              (without an over-approximation)
%
% Outputs:
%    pZsplit - cell array of split polyZonotopes
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
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: split, splitOneGen

% Author:       Niklas Kochdumper
% Written:      29-March-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% find all generators with a polynomial order of one
temp = sum(pZ.expMat,1);
ind = find(temp == 1);

% check if there are generators with a polynomial order of one that can be
% used for splitting
if ~isempty(ind)
    % determine the longest generator with polynomial order 1
    temp = sum(pZ.G(:,ind).^2,1);
    [~,indNew] = sort(temp,'descend');
    index = ind(indNew(1));
else
    % determine the longest generator
    temp = sum(pZ.G.^2,1);
    [~,indNew] = sort(temp,'descend');
    index = indNew(1);
end

% split the zonotope at the determined generator
if nargin == 2
   pZsplit = splitOneGen(pZ,index,varargin{1}); 
else
   pZsplit = splitOneGen(pZ,index); 
end

%------------- END OF CODE --------------