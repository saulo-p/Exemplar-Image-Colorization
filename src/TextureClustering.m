function clusters = TextureClustering(src_txt, lab_img, nClusters, PLOT)  
%% K-means clustering

if (isempty(lab_img))
  ab = [];
else
  ab = double(lab_img);
  nrows = size(ab,1);
  ncols = size(ab,2);

  ab = reshape(ab(:,:,2:3), nrows*ncols, 2); 
end
% repeat the clustering process N times to avoid local minima
[cluster_idx, C, ~, D] = kmeans([src_txt' ab], nClusters, 'Distance', 'sqEuclidean', ...
                          'Replicates', 5);

   
%% Cluster statistics

minD = min(D');

stds = zeros(1, nClusters);
sizes = zeros(1, nClusters);

for i = 1:nClusters
    idx = cluster_idx == i;
    
    stds(i) = std(minD(idx));
    sizes(i) = sum(idx);
end

%% Misclassification costs matrix

% mc_cost = squareform(pdist(C(:,1:2)));
% if (true)
%   for i = 1:size(mc_cost,2)
%     mc_cost(:,i) = mc_cost(:,i)/sum(mc_cost(:,i));
%   end
% else
%   mc_cost = mc_cost*((nClusters - 1)/norm(mc_cost));
% end

%% Output

clusters.idxs = cluster_idx;
clusters.centroids = C;
clusters.stds = stds;
clusters.cardin = sizes;
% clusters.mcCost = mc_cost;

end