function [] = StructHA2file(Data,functionName,resultpath)
% Function creates a new Matlab file containing a parallel hybrid
% automaton in CORA formate
%
%INPUT:
%    Data: Automaton in structHA format
%    functionName (optional): desired name of generated MATLAB function
%               (default: filename of source SX file, contained in Data)
%    resultpath (optional): target directory of generated files
%               (default: <transformer directory>/coramodels)
%OUTPUT:
%    void

if nargin >= 1
    % if a function name is given, overwrite the default name of Data
    Data.name = functionName;
end

if nargin < 3
    % if no resultpath is given, use "cora/models/SpaceExConverted"
    resultpath = strcat(coraroot,'/models/SpaceExConverted/',Data.name);
end

if ~exist(resultpath,'dir')
    mkdir(resultpath);
end
addpath(resultpath);

[Filename,Str] = data2ParallelHA(Data,resultpath);

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
    error('could not open output file %s',fname);
end

fwrite(fileID, Str);
fclose(fileID);

% ensure matlab detects new function
rehash path;

end

