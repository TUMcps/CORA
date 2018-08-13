function plotAsText(varargin)
% plotAsText - Plots vertices of a 2-dim projection as text
%
% Syntax:  
%    plotAsText(Rcont_zono,projectedDimensions,fileName);
%
% Inputs:
%    Rcont_zono - cell array of zonotope objects
%    projectedDimensions - dimensions that should be projected (optional) 
%    fileName - file name of the text file
%
% Outputs:
%    none
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      29-January-2014
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%If only one argument is passed
if nargin == 2
    Z = varargin{1};
    Rcont_zono = varargin{2};
    dimensions = [1,2];
    fileName = 'vertices';
    
%If two arguments are passed    
elseif nargin == 3
    Z = varargin{1};
    Rcont_zono = varargin{2};
    dimensions = varargin{3};
    fileName = 'vertices';
    
%If too many arguments are passed
elseif nargin == 4
    Z = varargin{1};
    Rcont_zono = varargin{2};
    dimensions = varargin{3};   
    fileName = varargin{4};
end

%open file
fid = fopen([fileName,'.txt'],'w');

for i = 1:length(Rcont_zono)
    %delete zero generators
    Z=deleteZeros(Rcont_zono{i});

    %Compute potential vertices
    V=vertices(Z,dimensions);
    V_mat = get(V,'V');
    
    %fprintf(fid, '%f\n', V_mat);
    %fprintf(fid, '%s\n', num2str(V_mat,10));
    writeMatrix(V_mat,fid)
end

%close file
fclose(fid);


function writeMatrix(M,fid)

%write each row
for iRow=1:(length(M(:,1)))
    if (length(M(1,:))-1)>0
        for iCol=1:(length(M(1,:))-1)
            str=num2str(M(iRow,iCol),10);
            str=[str,', '];
            %write in file
            fprintf(fid, '%s', str);
        end
    else
        iCol = 0; %for vectors
    end
    if iRow<length(M(:,1))
        %write last element
        str=num2str(M(iRow,iCol+1),10);
        %write in file
        fprintf(fid, '%s\n', str);
    else
        %write last element
        str=num2str(M(iRow,iCol+1),10);
        %write in file
        fprintf(fid, '%s\n\n', str);   
    end
end


%------------- END OF CODE --------------