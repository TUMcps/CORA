function [id,msg] = errConstructor(varargin)
% errConstructor - standardized error message format if the constructor
%    of an object receives wrong inputs, e.g., dimension mismatch
%    between inputs or wrong format; this functions accounts for all
%    mistakes which cannot be easily found by property validation
%
% Syntax:  
%    [id,msg] = errConstructor(varargin)
%
% Inputs:
%    infomsg - message containing information about error
%
% Outputs:
%    msg - error message
%    id - error identifier
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: noops

% Author:        Mark Wetzlinger
% Written:       19-March-2021
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

infomsg = [];
if ~isempty(varargin) && ischar(varargin{1})
    infomsg = varargin{1};
end

st = dbstack;
filename = st(2).file; % calling function name: constructor
classname = filename(1:end-2);

id = 'CORA:wrongInputInConstructor';
msg = ['Wrong input arguments for constructor of class ' classname '!'];
if ~isempty(infomsg)
    msg = [msg '\n   Information: ' infomsg]; 
end


%------------- END OF CODE --------------
