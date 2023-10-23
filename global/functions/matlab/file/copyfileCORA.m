function copyfileCORA(source,target)
% copyfileCORA - copies a file from a source to a destination; this
% function has been created since when using the built-in MATLAB function
% 'copyfile', the created file is not recognized by all MATLAB versions
% even when using addpath(genpath(path));
%
% Syntax:
%    copyfileCORA(source,target)
%
% Inputs:
%    source - source of the file
%    target - target of the file
%    belong
%
% Outputs:
%    -
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       12-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% open source file
fid_source = fopen(source,'r');

% create new file
fid_target = fopen(target,'w');

% scan content of source file line by line
while ~feof(fid_source)
    % read next line
    line = fgetl(fid_source);
    % write line into file
    fprintf(fid_target, '%s \n', line);
end

%close files
fclose(fid_source);
fclose(fid_target);

% ------------------------------ END OF CODE ------------------------------
