function xhits = scangenes(genes,numx)
%SCANGENES Scan a single multigene individual for all input variables and return a frequency vector.
%
%   XHITS = SCANGENES(GENES,NUMX) scans the cell array GENES for NUMX
%   input variables and returns the frequency vector XHITS.
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also GPMODELVARS

numgenes = length(genes);
xhits = zeros(1,numx);

%loop through inputs
for i=1:numx
    k1 = [];
    k2 = [];
    istr = ['x' int2str(i)];
    
    %loop through genes, look for current input
    for j=1:numgenes
        
        k1 = strfind(genes{j}, [istr ',']);
        k2 = strfind(genes{j},[istr ')']);
        
        %workaround for special case (trees containing a single terminal node)
        k3 = strfind(genes{j},istr);
        if ~isempty(k3) && numel(k3) == 1
            if length(istr) == length(genes{j});
                xhits(i) = xhits(i) + 1;
            end
        end
        
        if ~isempty(k1)
            xhits(i) = xhits(i) + length(k1);
        end
        
        if ~isempty(k2)
            xhits(i) = xhits(i) + length(k2);
        end
        
    end
    
end