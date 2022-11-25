function [ test ] = isFinalLocation(loc, finalLoc )
%ISFINALLOCATION Summary of this function goes here
%   Detailed explanation goes here

if(finalLoc == 0)
    test = false;
    
else
    [nFinalLocations, ~] = size(finalLoc);
    res = zeros(nFinalLocations,1);
    for iFinalLocations = 1:1:nFinalLocations
        
        res(iFinalLocations) = sum((loc == finalLoc(iFinalLocations,:)));
    end
    
    test = (sum(res)>0);
end

end

