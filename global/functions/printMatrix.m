function printMatrix(M,varargin)
% printMatrix - prints an matrix such that if one executes this command
%    in the workspace, this matrix would be created
%
% Syntax:  
%    printMatrix(M)
%
% Inputs:
%    M - matrix
%    accuracy - (optional) floating-point precision
%
% Outputs:
%    -
%
% Example: 
%    M = [2 3; -2 1];
%    printMatrix(M)

% Author:       Matthias Althoff
% Written:      01-November-2017
% Last update:  27-June-2018
%               04-Jan-2021
%               17-June-2022 (MW, parsing of accuracy)
% Last revision:---

%------------- BEGIN CODE --------------

% determine accuracy
% default accuracy
accuracy = '%4.3f%s';
if nargin == 2
    % numerical value
    if isnumeric(varargin{1})
        accuracy = varargin{1};
    % category
    elseif ischar(varargin{1})
        if strcmp(varargin{1},'high')
            accuracy = '%16.16f%s';
        end
    end
elseif nargin > 2
    throw(CORAerror('CORA:tooManyInputArgs',2));
end

% write first element
fprintf('%s\n','[ ...');

%write each row
for iRow=1:(length(M(:,1)))
    if (length(M(1,:))-1)>0
        for iCol=1:(length(M(1,:))-1)
            %write in workspace
            fprintf(accuracy, M(iRow,iCol), ', ');
        end
    else
        iCol = 0; %for vectors
    end
    if iRow<length(M(:,1))
        %write in workspace
        fprintf([accuracy,'\n'], M(iRow,iCol+1), '; ...');
    else
        %write in workspace
        fprintf([accuracy,'\n\n'], M(iRow,iCol+1), '];');   
    end
end

%------------- END OF CODE --------------
