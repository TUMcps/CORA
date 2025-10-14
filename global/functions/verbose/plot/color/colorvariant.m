function [color,alpha] = colorvariant(variant,color,varargin)
% colorvariant - provides 'light' and 'dark' color variants
%
% Syntax:
%    color = colorvariant(variant,color)
%    color = colorvariant(variant,color,alpha)
%
% Inputs:
%    variant - str, 'light', 'dark', 'none'
%    color - numeric, color triplet
%    alpha - numeric, factor for lightening/darkening. Default: 0.2
%
% Outputs:
%    color - numeric, updated color triplet
%    alpha - numeric, factor for lightening/darkening

% Authors:       Lukas Koller, Tobias Ladner
% Written:       11-September-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(2,3)
alpha = setDefaultValues({0.2},varargin);
inputArgsCheck({ ...
    {variant,'str',{'light','dark','none'}}; ...
    {color,'att','numeric'}; ...
    {alpha,'att','numeric'};
})

% Apply a color modification based on the postfix.
switch variant
    case 'light'
        % Combine the color with white.
        color = 1.*(1-alpha) + color.*alpha;
    case 'dark'
        % Combine the color with black.
        color = 0.*(1-alpha) + color.*alpha;
    case 'none'
        % Apply no modification.
    otherwise
end

end

% ------------------------------ END OF CODE ------------------------------
