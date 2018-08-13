function outputTest = test_parallelHybridAutomaton()

HA = testFlatHybridAutomaton();

PHA = testParallelHybridAutomaton();

rHA = get(HA,'trajectory');

rPHA = get(PHA,'trajectory');

locHA = rHA.loc;
locPHA = rPHA.loc;

dim2 = size(locHA,2);
dim1 = size(locPHA,1);

minDim = min(dim1,dim2);

for iCount = 1:1:minDim
    
    if (isequal(locPHA(iCount,:),[1,1]))
        res(iCount) = 1;
    elseif (isequal(locPHA(iCount,:),[1,2]))
        res(iCount) = 2;
    elseif (isequal(locPHA(iCount,:),[2,1]))
        res(iCount) = 4;
    elseif (isequal(locPHA(iCount,:),[2,2]))
        res(iCount) = 3;
    else
        error('test_parallelHybridAutomaton error');
    end

end

try
sumError = sum(locHA ~= res);
catch
sumError = 1;
end

outputTest = (sumError == 0);

end



