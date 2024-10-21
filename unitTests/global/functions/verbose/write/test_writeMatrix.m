function res = test_writeMatrix
% test_writeMatrix - unit tests for writing the contents of nD symbolic 
%    matrices to standalone files
%
% Syntax:
%    res = test_writeMatrix
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       13-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set path to this folder
path = mfilename('fullpath');
idxFilesep = strfind(path,filesep);
path = path(1:idxFilesep(end));
fname = 'test_writeMatrix_generatedfile';
fullname = [path filesep fname '.m'];

try

    % 2D
    M = sym([1 2 3; 4 5 6]);
    fid = fopen(fullname,'w');
    writeMatrix(fid,M,'M');
    fclose(fid);
    delete(fullname);
    
    % 3D
    M = sym(zeros(2,3,2));
    M(:,:,1) = sym([1 2 3; 4 5 6]);
    M(:,:,2) = sym([7 8 9; 10 11 12]);
    fid = fopen([path filesep fname '.m'],'w');
    writeMatrix(fid,M,'M');
    fclose(fid);
    delete(fullname);
    
    % 4D
    M = sym(reshape(1:120,[4,5,3,2]));
    fid = fopen([path filesep fname '.m'],'w');
    writeMatrix(fid,M,'M','Sparse',false);
    fclose(fid);
    delete(fullname);

catch ME
    fclose(fid);
    rethrow(ME);
end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
