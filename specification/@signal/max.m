function val = max(sig, int)
% max - get the maximum value within the given interval
%
% Syntax:
%    val = max(sig, int)
%
% Inputs:
%    sig - input signal
%    int - interval for the maximum value
%
% Outputs:
%    val - maximum value
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Benedikt Seidl
% Written:       23-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

arguments
    sig signal
    int interval
end

val = 0;

% find index where interval begins
for i=1:length(sig)
    if sig.time(i) >= infimum(int)
        val = max(val, sig.value(i));
        break
    end
end

% continue while interval lasts
for j=i:length(sig)-1
    if sig.time(j) <= supremum(int)
        val = max(val, sig.value(j+1));
    end
end

end

% ------------------------------ END OF CODE ------------------------------
