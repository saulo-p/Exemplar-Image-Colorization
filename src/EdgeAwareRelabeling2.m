function relabels = EdgeAwareRelabeling2(clusters, labels, labels_scores)
  relabels = zeros(size(labels));
  for ci = 1:length(clusters)
    cluster_labels = labels(clusters{ci});
    cluster_labels = cluster_labels(cluster_labels ~= -1);
    if (~isempty(cluster_labels))
      relabels(clusters{ci}) = mode(cluster_labels);
    else
      relabels(clusters{ci}) = -1;
    end
  end
  
  %Only changes if the original score of the new class is high enough;
  scores_idxs = sub2ind(size(labels_scores), (1:size(labels_scores,1))', relabels);
  relabels_scores = labels_scores(scores_idxs);
  
  %New label needs at least probability of 1/nClasses.
  confidence = relabels_scores > 1/size(labels_scores,2);
  relabels = relabels.*confidence + labels.*~confidence;
end

