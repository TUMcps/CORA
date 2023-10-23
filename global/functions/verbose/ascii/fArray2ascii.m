function fArray2ascii
% fArray2ascii - ???
%
% Syntax:
%    fArray2ascii
%
% Inputs:
%    - 
%
% Outputs:
%    - 
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Matthias Althoff
% Written:       23-May-2008
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%open .mat file
[FileName,PathName] = uigetfile();
cd(PathName);
file=load(FileName);
fArray=file.fArray;

%create filename
fName=strrep(FileName,'.mat','.txt');
%open file for writing
fid = fopen(fName,'w');

%get nr of angle, x and y segments
nrOfAngles=length(fArray.val);
nrOfXseg=length(fArray.val{1}(:,1));
nrOfYseg=length(fArray.val{1}(1,:));

%write to file
%write number of angles, x- and y-values first
fprintf(fid, '%4i %4i %4i \n', [nrOfAngles,nrOfXseg,nrOfYseg]);
%write segment lenght for angle, x-pos, y-pos
fprintf(fid, '%6i %6i %6i \n', [fArray.segLength.angle, fArray.segLength.x, fArray.segLength.y]);
%write values
for iAngle=1:nrOfAngles
    for iX=1:nrOfXseg
        for iY=1:nrOfYseg
            %save nr of rows, columns, nonzero elements
            fprintf(fid, '%4i %4i %4i %6i \n', [iAngle, iX, iY, fArray.val{iAngle}(iX,iY)]);
        end
    end
end
    
%close file
status = fclose(fid)

% ------------------------------ END OF CODE ------------------------------
