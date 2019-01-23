%% -----------------------------------------------------------------
% Beyond Lanscapes Exemplar-based Method
% Author: Saulo Pereira
%-------------------------------------------------------------------

input_file = 'single';

%% Input parameters
[IP, FP, OO] = InputAlgorithmParameters(input_file);
figs = GenerateFiguresList;

%% Input data (source and target images)
[source.image, target.image] = LoadImages(IP.sourceFile, IP.targetFile, IP.dataFolder);

%% Color space conversion
source.lab = rgb2lab(source.image);

if (OO.PLOT)
  ShowColorDistribution(source.image, source.lab);
end

%% Luminance Remapping
tgt_lab = rgb2lab(cat(3, target.image, target.image, target.image));
target.luminance = tgt_lab(:,:,1)/100;

source.luminance = luminance_remap(source.lab, target.luminance, IP.sourceFile == IP.targetFile);

%% Color Clustering (Automatic Labeling)
% Performs the clustering for classification.
if (IP.COLOR_CLUSTERING)
  disp('Source color clustering'); tic;
  clusters = ColorClustering(source.lab, IP.nClusters, IP.CL_CHANNELS, OO.PLOT);

  toc;
end

%%
[samples.idxs, samples.ab] = FullSampling(source.lab);
samples.sourceSize = size(source.luminance);

%% Feature extraction
disp('Feature extraction'); tic;

[target.fv, target.fvl] = FeatureExtraction(target.luminance, FP);
[samples.fv, samples.fvl] = FeatureExtraction(source.luminance, FP);
toc;

%% Superpixel extraction
disp('Superpixel extraction'); tic;
[source.sp, source.lin_sp, source.sp_centroids, source.nSuperpixels] = ...
  SuperpixelExtraction(source.luminance, IP.nSuperpixels, 'turbo');
[target.sp, target.lin_sp, target.sp_centroids, target.nSuperpixels] = ...
  SuperpixelExtraction(target.image, IP.nSuperpixels, 'turbo');
toc;

if (OO.PLOT)
  %Show superpixels
  figure(figs.TargetSP); imshow(imoverlay(target.image, boundarymask(target.sp, 4), 'w')); 
  title('Target superpixels');
  figure(figs.SourceSP); imshow(imoverlay(source.image, boundarymask(source.sp, 4), 'w')); 
  title('Source superpixels');
end

%>Superpixel Labeling
if (IP.COLOR_CLUSTERING)
  disp('Superpixel labeling'); tic;
  [source.sp_clusters, source.sp_chrom] = SuperpixelLabeling(source.lab, clusters.idxs, source.lin_sp, ...
    IP.LBL_MAJOR, OO.PLOT, source.sp, size(source.luminance));
end

%> Superpixel Feature Aggregation
disp('Superpixel feature aggregation'); tic;
[target.fv_sp, source.fv_sp] = SuperpixelsFeatures(source, samples, target);

target = rmfield(target, 'fv');
samples = rmfield(samples, 'fv');
toc;

%> Saliency Feature Computation
[ssi1, ssi2] = SaliencyFeature(source.luminance, source.sp, source.nSuperpixels);
[tsi1, tsi2] = SaliencyFeature(target.luminance, target.sp, target.nSuperpixels);

if (OO.PLOT)
  figure;
  subplot(2,1,1);
  imshow([ssi1 ssi2],[]);
  subplot(2,1,2);
  imshow([tsi1 tsi2],[]);
  title('Saliency Feature');
end

%Find unique superpixels indexes and concatenate their saliency values
%onto the feature vector.
[sp_idxs, src_idxs] = unique(source.lin_sp);
[~, tgt_idxs] = unique(target.lin_sp);
source.fv_sp = [source.fv_sp; ssi1(src_idxs)'; ssi2(src_idxs)'];
target.fv_sp = [target.fv_sp; tsi1(tgt_idxs)'; tsi2(tgt_idxs)'];
clear ssi1 ssi2 tti1 tti2;

%% Classification
disp('Classification in Feature Space'); tic;
TWO_STAGE = true; %TODO: colocar na lista de inputs

if (TWO_STAGE)
  %% Source -> Target: initial label
  [~, nb_s1_idxs, nb_s1_dists] = CombinedPDists(source.fv_sp, target.fv_sp, FP.feats1p);
  nb_s1_classes = source.sp_clusters(nb_s1_idxs);
  
  %Stage 1 classification 
  [s1_labels, ~, ~, ~, ~] = ...
    PredictSuperpixelsClassesKNN(nb_s1_classes, nb_s1_dists, IP.Kfs, IP.nClusters, clusters.mcCost);
  
  if (OO.PLOT)
    surf_labels_img = CreateLabeledImage(s1_labels, target.sp, size(target.image));

    figure; imshow(surf_labels_img, []); colormap jet;
    title('Pass-1 labels');
    drawnow;
    clear surf_labels_img
  end
  
  clear nb_s1_classes nb_s1_dists nb_s1_idxs;
  
  %% Target -> Target: label refinement
  [~, nb_s2_idxs, nb_s2_dists] = CombinedPDists(target.fv_sp, target.fv_sp, FP.feats2p);
  nb_s2_classes = s1_labels(nb_s2_idxs);

  %Intra-image (remove the smallest distance which is the element itself).
  nb_s2_classes = [nb_s2_classes(:,2:end) nb_s2_classes(:,1)];
%   nb_s2_idxs = [nb_s2_idxs(:,2:end) nb_s2_idxs(:,1)];
  nb_s2_dists = [nb_s2_dists(:,2:end) nb_s2_dists(:,1)];
  
  %s2_labels: kNNE with cost term.
    [~, s2_labels, ~, final_scores, ~] = ...
    PredictSuperpixelsClassesKNN(nb_s2_classes, nb_s2_dists, IP.Kfs, IP.nClusters, clusters.mcCost);
  
  if (OO.PLOT)
    two_stage_labels_img = CreateLabeledImage(s2_labels, target.sp, size(target.image));

    figure; imshow(two_stage_labels_img, []); colormap jet;
    title('Two-stage labels (cost kNNE)');
    drawnow;
  end
  
  clear nb_s2_classes nb_s2_dists nb_s2_idxs;
else
  %>Single-pass
  [~, nb_s0_idxs, nb_s0_dists] = CombinedPDists(source.fv_sp, target.fv_sp, ...
    FP.feats1p + FP.feats2p);
  nb_s_classes = source.sp_clusters(nb_s0_idxs);
  %kNN Equality
  [~, sing_labels, ~, ~, ~] = ...
    PredictSuperpixelsClassesKNN(nb_s_classes, nb_s0_dists, IP.Kfs, IP.nClusters, clusters.mcCost);

  if (OO.PLOT)
    sing_stage_labels_img = CreateLabeledImage(sing_labels, target.sp, size(target.image));

    figure; imshow(sing_stage_labels_img, []); colormap jet;
    title('Two-stage labels (cost kNNE)');
    drawnow;
  end
  
  clear nb_s0_idxs nb_s0_dists nb_s0_classes 
end


%>Color neighborhood
[~, neighbor_idxs, neighbor_dists] = CombinedPDists(source.fv_sp, target.fv_sp, FP.feats1p + FP.feats2p);
neighbor_classes = source.sp_clusters(neighbor_idxs);

%% Edge-Aware Labeling/Relabeling
disp('Edge-Aware Labeling/Relabeling'); tic;

%Clustering superpixels by adapted connected component labeling
[eaClusters, eaClustersImg] = EdgeAwareClustering(target);
%Relabeling of superpixels
relabeled = EdgeAwareRelabeling(eaClusters, s2_labels, final_scores);

if (OO.PLOT)
  relabeled_img = CreateLabeledImage(relabeled, target.sp, size(target.image));

  figure; imshow(relabeled_img, []); colormap jet;
  title('Relabels> [costkNNE]');
  drawnow;
end

toc;

%% Color transfer:
disp('Color transfer + File Save'); tic

[tgt_scribbled, scribbles_mask] = CopyClosestSuperpixelFromClassAvgScribble(source, target, ...
      neighbor_idxs, neighbor_classes, relabeled, IP.Kfs);
target.rgb = ColorPropagationLevin(lab2rgb(tgt_scribbled), target.luminance, scribbles_mask);
target.rgb_sat = IncreaseSaturation(target.rgb, 0.2);
imwrite(target.rgb_sat, ['./../results/' IP.sourceFile(1:end-2-3) '_final.png'], 'png');

toc;
