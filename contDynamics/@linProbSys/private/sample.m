function [obj] = sample(obj)
% sample - computes a sample of the system matrix
%
% Syntax:  
%    [obj] = sample(obj)
%
% Inputs:
%    obj - linIntSys object
%
% Outputs:
%    obj - linIntSys object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      16-May-2007 
% Last update:  15-June-2016
%               25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------

%get left and right system matrix
[Aleft,Aright]=interval(obj.A);

%convert to interval matrix
intervals=[reshape(Aleft,[],1),reshape(Aright,[],1)];

%build interval
I=interval(intervals(:,1), intervals(:,2));

%get vertices of interval
V = vertices(I);
W = unique(V', 'rows')';

%reshape to different As
[rows,cols]=size(Aleft);
for i=1:length(W(1,:))
    obj.sample.A{i}=reshape(W(:,i),rows,cols);
end
    
%------------- END OF CODE --------------