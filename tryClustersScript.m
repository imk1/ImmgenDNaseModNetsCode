data = importdata('nonTFExpressionsFiltered5');
save nonTFExpressionsFiltered5 data
load('cellTypesTrain.mat')
nonTFDataTrain = data.data(:,cellTypesTrain);
save nonTFDataTrain5 nonTFDataTrain
load('cellTypesTest.mat')
nonTFDataTest = data.data(:,cellTypesTest);
save nonTFDataTest5 nonTFDataTest

[L, C] = kmeansPlusPlus(nonTFDataTrain', 100);
idx = kmeans(nonTFDataTrain, 100, 'start', C');
sils = silhouette(nonTFDataTrain, idx);
mean(sils)
median(sils)
clustSils = zeros(max(idx), 1);
for i = 1:max(idx)
    clustLocs = find(idx == i);
    currentSils = sils(clustLocs);
    clustSils(i) = mean(currentSils);
end
length(find(clustSils > .5))

[L, C] = kmeansPlusPlus(nonTFDataTrain', 50);
idx = kmeans(nonTFDataTrain, 50, 'start', C');
sils = silhouette(nonTFDataTrain, idx);
mean(sils)
median(sils)
clustSils = zeros(max(idx), 1);
for i = 1:max(idx)
    clustLocs = find(idx == i);
    currentSils = sils(clustLocs);
    clustSils(i) = mean(currentSils);
end
length(find(clustSils > .5))

[L, C] = kmeansPlusPlus(nonTFDataTrain', 10);
idx = kmeans(nonTFDataTrain, 10, 'start', C');
sils = silhouette(nonTFDataTrain, idx);
mean(sils)
median(sils)
clustSils = zeros(max(idx), 1);
for i = 1:max(idx)
    clustLocs = find(idx == i);
    currentSils = sils(clustLocs);
    clustSils(i) = mean(currentSils);
end
length(find(clustSils > .5))

nonTFDataCorr = corrcoef(nonTFDataTrain);
mean(mean(nonTFDataCorr))
min(min(nonTFDataCorr))
count = 0;
for i = 1:147
    for j = i+1:147
        if nonTFDataCorr(i,j) < .5
            count = count + 1;
        end
    end
end
total = nchoosek(147, 2);
count/total

idx = hierarchicalClusterDivisive(nonTFDataTrain, 10, 10);
sils = silhouette(nonTFDataTrain, idx);
mean(sils)
median(sils)
clustSils = zeros(max(idx), 1);
for i = 1:max(idx)
    clustLocs = find(idx == i);
    currentSils = sils(clustLocs);
    clustSils(i) = mean(currentSils);
end
length(find(clustSils > .5))

idx = hierarchicalClusterDivisive(nonTFDataTrain, 100, 10);
sils = silhouette(nonTFDataTrain, idx);
mean(sils)
median(sils)
clustSils = zeros(max(idx), 1);
for i = 1:max(idx)
    clustLocs = find(idx == i);
    currentSils = sils(clustLocs);
    clustSils(i) = mean(currentSils);
end
length(find(clustSils > .5))

idx = hierarchicalClusterDivisive(nonTFDataTrain, 10, 100);
sils = silhouette(nonTFDataTrain, idx);
mean(sils)
median(sils)
clustSils = zeros(max(idx), 1);
for i = 1:max(idx)
    clustLocs = find(idx == i);
    currentSils = sils(clustLocs);
    clustSils(i) = mean(currentSils);
end
length(find(clustSils > .5))

idx = hierarchicalClusterDivisive(nonTFDataTrain, 100, 100);
mean(sils)
median(sils)
clustSils = zeros(max(idx), 1);
for i = 1:max(idx)
    clustLocs = find(idx == i);
    currentSils = sils(clustLocs);
    clustSils(i) = mean(currentSils);
end
length(find(clustSils > .5))

idx = clusterdata(nonTFDataTrain, 100);
sils = silhouette(nonTFDataTrain, idx);
mean(sils)
median(sils)
clustSils = zeros(max(idx), 1);
for i = 1:max(idx)
    clustLocs = find(idx == i);
    currentSils = sils(clustLocs);
    clustSils(i) = mean(currentSils);
end
length(find(clustSils > .5))

idx = clusterdata(nonTFDataTrain, 10);
sils = silhouette(nonTFDataTrain, idx);
mean(sils)
median(sils)
clustSils = zeros(max(idx), 1);
for i = 1:max(idx)
    clustLocs = find(idx == i);
    currentSils = sils(clustLocs);
    clustSils(i) = mean(currentSils);
end
length(find(clustSils > .5))

G = makeAssociationKernel(nonTFDataTrain', 100);
[idx,netsim,dpsim,expref]=apcluster(G,0);
indexes = unique(idx);
idxPlus = zeros(length(idx), 1);
for i = 1:length(idx)
    index = find(indexes == idx(i));
    idxPlus(i) = index;
end
sils = silhouette(nonTFDataTrain, idxPlus);
mean(sils)
median(sils)
clustSils = zeros(max(idxPlus), 1);
for i = 1:max(idxPlus)
    clustLocs = find(idxPlus == i);
    currentSils = sils(clustLocs);
    clustSils(i) = mean(currentSils);
end
length(find(clustSils > .5))

[coeff,score,latent] = princomp(nonTFDataTrain);
var = cumsum(latent)./sum(latent);
[L, C] = kmeansPlusPlus(score(:,1:16)', 20);
idx = kmeans(score(:,1:16), 20, 'start', C');
sils = silhouette(score(:,1:16), idx);
mean(sils)
median(sils)
clustSils = zeros(max(idx), 1);
for i = 1:max(idx)
    clustLocs = find(idx == i);
    currentSils = sils(clustLocs);
    clustSils(i) = mean(currentSils);
end
length(find(clustSils > .5))

idxArray = zeros(size(nonTFDataTrain));
for i = 1:size(nonTFDataTrain, 1)
    idx = kmeans(nonTFDataTrain(i, :), 2);
    idxArray(i, :) = idx';
end
sils = zeros(size(nonTFDataTrain));
for i = 1:size(nonTFDataTrain, 1)
    sils(i, :) = silhouette(nonTFDataTrain(i, :)', idxArray(i, :)')';
end
mean(mean(sils))
clusterableGenes = [];
for i = 1:size(sils,1)
    if mean(sils(i, :)) > .5
        clusterableGenes = vertcat(clusterableGenes, i);
    end
end
for i = 1:size(idxArray, 1)
    if mean(nonTFDataTrain(find(idxArray(i,:) == 1))) > mean(nonTFDataTrain(find(idxArray(i,:) == 2)))
        for j = 1:size(idxArray, 2)
            if idxArray(i, j) == 1
                idxArray(i, j) = 2;
            else
                idxArray(i, j) = 1;
            end
        end
    end
end
[L, C] = kmeansPlusPlus(idxArray(clusterableGenes, :)', 20);
idx = kmeans(idxArray(clusterableGenes, :), 20, 'start', C');
silsPlus = silhouette(idxArray(clusterableGenes, :), idx);
mean(silsPlus)
median(silsPlus)
clustSils = zeros(20, 1);
for i = 1:20
    clustLocs = find(idx == i);
    currentSils = silsPlus(clustLocs);
    clustSils(i) = mean(currentSils);
end
length(find(clustSils > .5))

nonTFTrainDataNorm = zeros(size(nonTFDataTrain));
for i = 1:size(nonTFDataTrain, 1)
    nonTFDataTrainNorm(i, :) = (nonTFDataTrain(i,:) - mean(nonTFDataTrain(i, :))) / std(nonTFDataTrain(i,:));
end
[L, C] = kmeansPlusPlus(nonTFDataTrainNorm', 20);
idx = kmeans(nonTFDataTrainNorm, 20, 'start', C');
sils = silhouette(nonTFDataTrainNorm, idx);
mean(sils)
median(sils)
clustSils = zeros(20, 1);
for i = 1:20
    clustLocs = find(idx == i);
    currentSils = sils(clustLocs);
    clustSils(i) = mean(currentSils);
end
length(find(clustSils > .5))