function I = interval(obj)
% interval - encloses a mptPolytope by an interval
%
% Syntax:  
%    I = interval(obj)
%
% Inputs:
%    obj - mptPolytope object
%
% Outputs:
%    I - interval 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      01-February-2011
% Last update:  30-July-2016
% Last revision:---

%------------- BEGIN CODE --------------

%obtain bounding box in halfspace representation
try %MPT3
    B = outerApprox(obj.P);
    
    %get C matrix, d vector
    Hfull = B.H;
    C = Hfull(:,1:end-1);
    d = Hfull(:,end);
catch %MPT2
    B = bounding_box(obj.P);
    
    %get C matrix, d vector
    [C,d] = double(B);
end

nColumns = size(C,2);
leftLim = zeros(nColumns,1);
rightLim = zeros(nColumns,1);

for iColumn = 1:1:nColumns
    
    [rowInf,~]=find(C(:,iColumn)==-1);

    if(~isempty(rowInf))
        leftLim(iColumn) = -d(rowInf);
    else
        leftLim(iColumn) = -Inf;
    end
    
    [rowSup,~]=find(C(:,iColumn)==1);
    if(~isempty(rowSup))
        rightLim(iColumn) = d(rowSup);
    else
        rightLim(iColumn) = Inf;
    end
end


%instantiate interval
% check against numerical issues:
if all(withinTol(leftLim,rightLim))
    % choose one of the limits to avoid leftLim > rightLim jitter
    I = interval(leftLim);
else
    I = interval(leftLim,rightLim);
end




%------------- END OF CODE --------------