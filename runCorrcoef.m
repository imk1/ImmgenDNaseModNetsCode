function dataCorrMat = runCorrcoef(dataMatFileName, dataCorrMatFileName)
% Makes a correlationMatrix for data

dataMat = importdata(dataMatFileName);
dataCorrMat = corrcoef(dataMat);

save(dataCorrMatFileName, 'dataCorrMat');