function [ bindMatrice ] = matricesRep( obj )
%MATRICESREP Summary of this function goes here
%   Detailed explanation goes here


if(~isempty(obj))
    bindMatrice = zeros(1,3);
    
    bindMatrice(1,1) = obj.input;
    bindMatrice(1,2) = obj.component;
    bindMatrice(1,3) = obj.state;
else
    bindMatrice = [];
end

end

