function outputTest = test_SX2CORA()
%TEST_TRANSFORMER tests the transformer from SpaceEx to CORA
%   Detailed explanation goes here

% try
    obj = SX2CORATransformer();
    
    mFile = obj.generate('parallelHeaterExample');
    
    HA = testFlatHybridAutomaton();
    
    PHA = feval(mFile);
    
    testOptions;
    
    PHA = simulate(PHA,options);
    
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
            error('test1 error');
        end
        
    end
    
    try
        sumError = sum(locHA ~= res);
    catch
        warning('bug in simulation of the generated function');
        sumError = 1;
    end

    outputTest = (sumError == 0);
    
% catch
%     warning('bug in generation process of the transformer');
%     outputTest = 0;
% end

end

