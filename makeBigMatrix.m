%function EuclidDist = makeBigMatrix(dataFileName)
function EuclidDist = makeBigMatrix(S)
%data = importdata(dataFileName);
%S = sparse(data);
%clear data
%EuclidDist = [];
EuclidDist = zeros(size(S,1), size(S,1));
t1 = tic;
 for i = 1:size(S, 1)
    %t2 = tic;
    %display(['Processing ' num2str(i) '/' num2str(size(S, 1)) ])
    if mod(i, 500) == 0
        i
    end
    row = -sqrt(sum(((ones(size(S,1), 1)*S(i,:)) - S(:,:)).^2, 2));
%     row = zeros(1, size(S, 1));
%     for j = 1:size(S, 1)
%         dist = -norm(S(i,:) - S(j,:));
%         row(j) = dist;
%     end
    %EuclidDist = vertcat(EuclidDist, row');
    EuclidDist(i,:) = row;
    clear row
    %EuclidDist = sparse(EuclidDist);
    %toc(t2)
end
toc(t1)

save ../April2012Release/EucildDistNonTFExprNorm EuclidDist
%save ../PWMExpCutNormScores/EuclidDistMatmCD19+-DS13845 EuclidDist