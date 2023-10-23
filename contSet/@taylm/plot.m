function han = plot(tay,varargin)
% plot - plots an over-approximative projection of a Taylor model
%
% Syntax:
%    han = plot(tay)
%    han = plot(tay,dims)
%    han = plot(tay,dims,type)
%
% Inputs:
%    tay - taylm object
%    dims - (optional) dimensions for projection
%    type - (optional) plot settings (LineSpec and Name-Value pairs)
%           additional Name-Value pairs:
%               <'Splits',splits> - number of splits for refinement
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/plot

% Authors:       Tobias Ladner
% Written:       09-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1. convert to polyZonotope
pZ = polyZonotope(tay);

% 2. plot polyZonotope
han = plot(pZ,varargin{:});

% 3. clear han
if nargout == 0
    clear han;
end

end

% ------------------------------ END OF CODE ------------------------------
