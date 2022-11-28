function StructHA2file(Data,functionName,resultpath)
% StructHA2file - Function creates a new Matlab file containing a parallel
%    hybrid automaton
%
% Syntax:
%    StructHA2file(Data,functionName,resultpath)
%
% Inputs:
%    Data - Automaton in structHA format
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

% Author:       ???
% Written:      ???
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if nargin >= 1
    % if a function name is given, overwrite the default name of Data
    Data.name = functionName;
end

% Matlab does not allow '-' in file names, so we replace it with '_'
if contains(Data.name,'-')
    warning("Matlab does not allow '-' in file names, all occurrences of '-' are replaced with '_'!");
    Data.name = strrep(Data.name,'-','_');
end

if nargin < 3
    % if no resultpath is given, use "cora/models/SpaceExConverted"
    resultpath = [CORAROOT filesep 'models' filesep 'SpaceExConverted' filesep Data.name];
end

if ~exist(resultpath,'dir')
    mkdir(resultpath);
end
addpath(resultpath);

if length(Data.Components) == 1 && ...
        length(Data.Components{1,1}.States) == 1 && ...
        isempty(Data.Components{1,1}.States.Trans)
    % conversion to linear or non-linear system
    [Filename,Str] = data2NonLinSys(Data,resultpath);
else
    % conversion to hybrid system
    [Filename,Str] = data2ParallelHA(Data,resultpath);
end

fname = strcat(resultpath,'/',Filename,'.m');
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

if fileID<0
    throw(CORAerror('CORA:converterIssue',['Could not open output file ' fname]));
end

fwrite(fileID, Str);
fclose(fileID);

% ensure matlab detects new function
rehash path;

disp("----------------------StructHA2file COMPLETE---------------------");

%------------- END OF CODE --------------