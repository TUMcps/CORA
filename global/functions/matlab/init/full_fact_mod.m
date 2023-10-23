function des_mat = full_fact_mod(levels)
% full_fact_mod - gives full factorial design matrix for any levels  
%    more than 2 of any number of variables (minimum 2)
%
% Syntax:
%    des_mat = full_fact_mod(x1,x2,...);
%
% Inputs:
%    x1, x2, ... - variable levels either in row or column vector.
%                 These are not number of levels, but the levels itself.
%
% Outputs:
%    des_mat - mxn
%           m = total number of designs (product of all levels) and  
%           n = number of variables  
%           The first column shows all the first variable levels,
%           the second column shows second variable levels and so on. 
%
% Example: 
%    x1=[-1 1]; x2=[100:100:300];
%    des_mat = full_fact(x1,x2)

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Bhaskar Dongare, Matthias Althoff
% Written:       17-November-2008
% Last update:   08-April-2016 (MA)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~all(levels > 1)
    throw(CORAerror('CORA:wrongValue','first',...
        'Each variable should have minimum 2 levels'));
end

% Total number of design points  
total=prod(levels);
%Initilization of output matrix
des_mat=zeros(0,1);
% Loop for full factorial points
for i=1:length(levels)
    if i~=1 && i~=length(levels)
        temp=zeros(0,1);
        for j=1:prod(levels(1:i-1))
            temp1=repmat(1:levels(i),prod(levels(i+1:end)),1);
            temp1=sortrows(temp1);
            temp=[temp; temp1];
        end
    elseif i==length(levels)
        temp=repmat((1:levels(i))',total/levels(i),1);
    else
       temp=repmat((1:levels(i))',total/levels(i),1);
       temp=sortrows(temp);
    end
    des_mat=[des_mat temp];
end

% ------------------------------ END OF CODE ------------------------------
