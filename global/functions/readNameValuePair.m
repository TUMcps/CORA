function [NVpairs,value] = readNameValuePair(NVpairs,name,varargin)
% readNameValuePair - searches through plotOptions to find name from name-value
%    pair, returns corresponding value and plotOptions where searched
%    name-value pair is deleted
%
% Syntax:  
%    [plotOptions,value] = readNameValuePair(plotOptions,name)
%    [plotOptions,value] = readNameValuePair(plotOptions,name,check)
%    [plotOptions,value] = readNameValuePair(plotOptions,name,check,def)
%
% Inputs:
%    NVpairs - Name-Value pairs
%    name - name of name-value pair
%    check - (optional) functions that value should be check with
%    def - (optional) default values
%
% Outputs:
%    NVpairs - LineSpecification options + Name-Value pairs
%    value - value of name-value pair
%
% Example: 
%    plotOptions = {'r','Filled',true};
%    [plotOptions,value] = ...
%        readNameValuePair(plotOptions,'Filled','islogical');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: global/functions/readPlotOptions, polyZonotope/plot

% Author:        Mark Wetzlinger, Niklas Kochdumper
% Written:       15-July-2020 
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% default empty value (name not found)
if nargin >= 4
    value = varargin{2};
else
    value = [];
end

% check every second entry
for i=1:2:length(NVpairs)-1
    
    % has to be a char
    if ischar(NVpairs{i})
        
        if strcmp(NVpairs{i},name)
            % name found
            value = NVpairs{i+1};
            
            % check whether name complies with check
            if nargin >= 3
                if ~feval(varargin{1},value)
                    error("Invalid assignment for " + name);
                end
            end
            
            % empty the corresponding cells
            NVpairs{i} = []; NVpairs{i+1} = [];
            % ...and delete empty cells
            NVpairs = NVpairs(~cellfun('isempty',NVpairs));
            break
        end 
    end 
end

%------------- END OF CODE --------------
