function [pZsplit,factor] = splitLongestGen(pZ,varargin)
% splitLongestGen - Splits the longest generator dependent generator with a 
%    polynomial order of 1 for a polynomial zonotope
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
%    factor - identifier of the dependent factor that is split
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%
%    temp = splitLongestGen(pZ);
%    pZsplit1 = splitLongestGen(temp{1});
%    pZsplit2 = splitLongestGen(temp{2});
%
%    plot(pZ,[1,2],'FaceColor','r');
%
%    figure; hold on;
%    plot(pZsplit1{1},[1,2],'FaceColor','b');
%    plot(pZsplit1{2},[1,2],'FaceColor','g');
%    plot(pZsplit2{1},[1,2],'FaceColor','m');
%    plot(pZsplit2{2},[1,2],'FaceColor','c');
%
% Other m-files required: reduce
% Subfunctions: none
% MAT-files required: none
%
% See also: split, splitDepFactor

% Authors:       Niklas Kochdumper
% Written:       29-March-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% determine longest generator
len = sum(pZ.G.^2,1);
[~,ind] = max(len);

% find factor with the largest exponent
[~,factor] = max(pZ.E(:,ind));
factor = pZ.id(factor);

% split the zonotope at the determined generator
pZsplit = splitDepFactor(pZ,factor, varargin{:}); 

% ------------------------------ END OF CODE ------------------------------
