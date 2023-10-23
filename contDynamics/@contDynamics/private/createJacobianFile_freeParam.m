function createJacobianFile_freeParam(Jdyn,path,name)
% createJacobianFile_freeParam - generates an mFile that allows to compute
%    the jacobian for a certain state, input, and parameter value
%
% Syntax:
%    createJacobianFile_freeParam(obj)
%
% Inputs:
%    Jdyn - jacobians
%    path - path where the function should be created
%    name - function name for the Jacobian computation
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
% See also: ---

% Authors:       Matthias Althoff
% Written:       07-June-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% open file
fid = fopen([path filesep name '.m'],'w');
try

fprintf(fid, '%s\n\n', ['function [A,B]=',name,'(x,u,p)']);

% SYSTEM MATRIX
% write "A=["
fprintf(fid, '%s', 'A=[');
% write rest of matrix
writeMatrix(Jdyn.x,fid);

% INPUT MATRIX
% write "B=["
fprintf(fid, '%s', 'B=[');
% write rest of matrix
writeMatrix(Jdyn.u,fid);

catch ME
    % close file
    fclose(fid);
    
    rethrow(ME);
end

% close file
fclose(fid);

% ------------------------------ END OF CODE ------------------------------
