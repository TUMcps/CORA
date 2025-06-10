function res = closePrint(fid,closefid)
% closePrint - closes opened file if necessary
%
% Syntax:
%    res = closePrint(fid,closefid)
%
% Inputs:
%    fid - numeric, file identifier
%    closefid - logical, whether to close fid
%
% Outputs:
%    res - logical
%
%
% See also: initPrint

% Authors:       Tobias Ladner
% Written:       20-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(2,2);
res = true;

% close fid if necessary
if closefid
    res = ~fclose(fid);
end

end

% ------------------------------ END OF CODE ------------------------------
