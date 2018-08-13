function [error, errorInt, errorInt_x, errorInt_y, R_y] = linError_mixed_specialTOP(obj, options, R, Verror_y)
% linError - computes the linearization error
%
% Syntax:  
%    [obj] = linError(obj,options)
%
% Inputs:
%    obj - nonlinear DAE system object
%    options - options struct
%    R - actual reachable set
%
% Outputs:
%    error - zonotope overapproximating the linearization error
%    errorInt - interval overapproximating the linearization error
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff
% Written:      21-November-2011
% Last update:  23-May-2013
% Last revision:---

%------------- BEGIN CODE --------------

nrOfSplits = 4;
splits = 0;
Rnew{1} = R;

%split direction
splitDir = [1; zeros(obj.dim - 1,1)];

while splits < nrOfSplits
    Rsplit = split(Rnew{splits+1}, splitDir);
    
    %remove current Rnew
    Rnew(splits+1) = [];
    
    %add split sets
    Rnew{splits+1} = Rsplit{1};
    Rnew{splits+2} = Rsplit{2};
    
    %increment number of splits
    splits = splits + 1;
end

for iSet = 1 : length(Rnew)
    % compute individual linearization error
    [error{iSet}, errorInt{iSet}, errorInt_x{iSet}, errorInt_y{iSet}, R_y{iSet}] = linError_mixed_noInt(obj, options, Rnew{iSet}, Verror_y);
end


%------------- END OF CODE --------------