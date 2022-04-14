function printMatrix(M, varargin)
% printMatrix - prints an matrix such that if one executes this command
% in the workspace, this matrix would be created
%
% Syntax:  
%    printMatrix(M)
%
% Inputs:
%    no
%
% Outputs:
%    -
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      01-November-2017
% Last update:  27-June-2018
%               04-Jan-2021
% Last revision:---


%------------- BEGIN CODE --------------

% determine accuracy
if nargin == 1
    accuracy = '%4.3f%s';
elseif nargin == 2
    % numerical value
    if isnumeric(varargin{1})
        accuracy = varargin{1};
    % category
    elseif ischar(varargin{1})
        if strcmp(varargin{1},'high')
            accuracy = '%16.16f%s';
        end
    else
        % default accuracy
        accuracy = '%4.3f%s';
    end
else
    % default accuracy
    accuracy = '%4.3f%s';
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

