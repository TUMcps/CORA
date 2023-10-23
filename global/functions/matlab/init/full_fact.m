function des_mat = full_fact(varargin)
% full_fact - gives full factorial design matrix for any levels  
%    more than 2 of any number of variables (minimum 2)
%
% Syntax:
%    des_mat = full_fact(x1,x2,x3);
%    des_mat = full_fact([-1 1],[100 200 300],[1:4]);
%
% Inputs:
%    x1,x2,x3 - variable levels either in row or column vector. These are
%               not number of levels but the levels itself
%
% Outputs:
%    des_mat - mxn matrix where
%       m = total number of designs (product of all levels) and  
%       n = number of variables  
%            The first column shows all the first variable levels, second 
%            column shows second variable levels and so on.
%
% Example:
%    x1=[-1 1];x2=[100:100:300];
%    des_mat = full_fact(x1,x2) % OR
%    des_mat = full_fact([-1 1],[100:100:300])
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Bhaskar Dongare
% Written:       17-November-2008
% Last update:   17-June-2022 (MW, formatting)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

for i=1:nargin
    levels(i)=length(varargin{i});
    % Converting row vector to column vector
    if size(varargin{i},1)==1
        varargin{i}=varargin{i}';
    end
end
% Check number of variables and levels of each variable
if nargin<2
    throw(CORAerror('CORA:notEnoughInputArgs',2));
end
if ~all(levels >1)
    throw(CORAerror('CORA:wrongValue','first',...
        'Each variable should have minimum 2 levels'));
end

% Total number of design points  
total=prod(levels);
%Initilization of output matrix
des_mat=zeros(0,1);
% Loop for full factorial points
for i=1:nargin
    if i~=1 && i~=nargin
        temp=zeros(0,1);
        for j=1:prod(levels(1:i-1))
            temp1=repmat(varargin{i},prod(levels(i+1:end)),1);
            temp1=sortrows(temp1);
            temp=[temp; temp1];
        end
    elseif i==nargin
        temp=repmat(varargin{i},total/levels(i),1);
    else
       temp=repmat(varargin{i},total/levels(i),1);
       temp=sortrows(temp);
    end
    des_mat=[des_mat temp];
end

% ------------------------------ END OF CODE ------------------------------
