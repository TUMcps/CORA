function han = plot(sig,varargin)
% plot - plots a signal over its time axis
%
% Syntax:
%    han = plot(sig)
%
% Inputs:
%    sig - signal object
%
% Outputs:
%    han - handle to the graphics object
%
% Example:
%    sig = signal([1.2 2.3 3.0], [true false true]);
%    plot(sig)
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Benedikt Seidl
% Written:       16-August-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% duplicate every entry in the array
time = repelem(sig.time, 2);
value = repelem(sig.value, 2);

han = plot([0 time(1:end-1)], value, varargin{:});

if nargout == 0
    clear han;
end

end

% ------------------------------ END OF CODE ------------------------------
