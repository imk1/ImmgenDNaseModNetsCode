function cobj = runClustergram(dataMatrixFileName, clustergramFileName)
% Creates a clustergram

DMobj = importdata(dataMatrixFileName);
cobj = clustergram(DMobj);

save(clustergramFileName, 'cobj');