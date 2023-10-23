function StructHA2file(data,varargin)
% StructHA2file - Function creates a new Matlab file containing a parallel
%    hybrid automaton
%
% Syntax:
%    StructHA2file(data,functionName,resultpath)
%
% Inputs:
%    data - Automaton in structHA format
%    functionName (optional) - desired name of generated MATLAB function
%                 (default: filename of source SX file, contained in Data)
%    resultpath (optional) - target directory of generated files
%               (default: <transformer directory>/coramodels)
%
% Outputs:
%    ---
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       ???
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check for too many input arguments
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% set default values
[filename,resultpath] = setDefaultValues({data.name,...
    [CORAROOT filesep 'models' filesep 'SpaceExConverted' filesep data.name]},...
    varargin);

% MATLAB does not allow '-' in file names, so we replace it with '_'
if contains(filename,'-')
    warning("Matlab does not allow '-' in file names, " + ...
        "all occurrences of '-' are replaced with '_'!");
    filename = strrep(filename,'-','_');
end

% directory for resulting model file
if ~exist(resultpath,'dir')
    mkdir(resultpath);
end
addpath(resultpath);

if length(data.components) == 1 && ...
        length(data.components{1,1}.States) == 1 && ...
        isempty(data.components{1,1}.States.Trans)
    % conversion to linear or non-linear continuous system
    str = data2NonLinSys(data,filename,resultpath);
else
    % conversion to hybrid system
    str = data2ParallelHA(data,filename,resultpath);
end

% open file in path
fname = strcat(resultpath,'/',filename,'.m');
fileID = fopen(fname, 'w');

%FOR DEBUG (prevents overwriting of existing files)
% if ~exist(fname,'file')
%     fileID = fopen(fname, 'w');
% else
%     %if filename is already used, append "_xx"
%     for i = 2:99
%         fname = sprintf('%s/%s_%02d.m',resultpath,Filename,i);
%         if ~exist(fname,'file')
%             fileID = fopen(fname, 'w');
%             break;
%         end
%     end
% end
%END DEBUG

% file could not be opened
if fileID<0
    throw(CORAerror('CORA:converterIssue',['Could not open output file ' fname]));
end

% write entire string defining the CORA model into the file
fwrite(fileID, str);
% close the file
fclose(fileID);

% ensure matlab detects new function
rehash path;

disp("----------------------StructHA2file COMPLETE---------------------");

% ------------------------------ END OF CODE ------------------------------
