function [res] = test_nnLinearLayer_nnLinearLayer()
% test_nnLinearLayer_nnLinearLayer - tests constructor of nnLinearLayer
%
% Syntax:  
%    res = test_nnLinearLayer_nnLinearLayer()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Tobias Ladner
% Written:      23-November-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% check simple example
W = rand(4,3); 
b = rand(4,1);
layer = nnLinearLayer(W, b);

name = "TestLayer";
layer = nnLinearLayer(W, b, name);
if layer.name ~= name
    res = false;
end

% wrong input
try
    layer = nnLinearLayer();
    % should have thrown an error
    res = false;
catch
end
try
    layer = nnLinearLayer(W);
    % should have thrown an error
    res = false;
catch
end
try
    layer = nnLinearLayer(b);
    % should have thrown an error
    res = false;
catch
end
try
    layer = nnLinearLayer('TestLayer');
    % should have thrown an error
    res = false;
catch
end

% dimension missmatch
W = rand(4,3); 
b = rand(10,1);
try
    layer = nnLinearLayer(W, b);
    % should have thrown an error
    res = false;
catch
end


end

%------------- END OF CODE --------------

