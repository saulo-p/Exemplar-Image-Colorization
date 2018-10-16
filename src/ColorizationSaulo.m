%% -----------------------------------------------------------------
% Exemplar-based colorization algorithm
% Author: Saulo Pereira
%-------------------------------------------------------------------

clearvars -except batch_out batch_folder input_folder input_path list in_i;
% FP = [0.5,0.5,0.1,0.25,0.25,2,1,0,0.1];
input_file = 'default';
if(~exist('batch_folder')); batch_folder = ''; end
if(~exist('batch_out')); batch_out = 'TESTS/'; end
  
%Structure to hold the outputs of the experiment
exp_labels = {};
exp_cols = {};

% outTypes = {'kNNm', 'kNNE', 'costkNNE', 'costkNNER', 'match'};
outTypes = {'p1', 'p2kNNm', 'p2kNNE', 'pSkNNE', 'pSmatch', 'p2kNNEcost', 'costkNNER'};
fieldValues = [];
for i = 1:numel(outTypes)
  fieldValues = [fieldValues char(39) outTypes{i} char(39) ',' num2str(i) ','];
end
call = ['struct(' fieldValues(1:end-1) ')'];
outLabels = eval(call);

%% Input parameters
[IP, FP, OO] = InputAlgorithmParameters(input_file);
figs = GenerateFiguresList;

%% Input data (source and target images)
[source.image, target.image] = LoadImages(IP.sourceFile, IP.targetFile, IP.dataFolder);

%% Color space conversion
source.lab = rgb2lab(source.image);

if (OO.PLOT && false)
  ShowColorDistribution(source.image, source.lab);
end

%% Luminance Remapping
tgt_lab = rgb2lab(cat(3, target.image, target.image, target.image));
target.luminance = tgt_lab(:,:,1)/100;

source.luminance = luminance_remap(source.lab, target.luminance, IP.sourceFile == IP.targetFile);

%% Color Clustering (Automatic Labeling)
% Performs the clustering for sampling and/or classification.
if (IP.COLOR_CLUSTERING)
  disp('Source color clustering'); tic;
  clusters = ColorClustering(source.lab, IP.nClusters, IP.CL_CHANNELS, OO.PLOT);

  toc;
end

%% Source sampling
disp('Source image sampling'); tic;

switch IP.SAMPLE_METHOD
  case 0
  %No sampling:
  [samples.idxs, samples.ab] = FullSampling(source.lab);

  case 1
  %Jittered sampling:
  [samples.idxs, samples_ab] = JitterSampleIndexes(source.lab, IP.nSamples);
  samples.idxs = [samples.idxs(2,:); samples.idxs(1,:)];
  samples.lin_idxs = sub2ind(size(source.luminance), samples.idxs(1,:), samples.idxs(2,:))';
  samples.ab = samples_ab(2:3,:);

  case 2
  %Clustered sampling:
  samples = ClusteredSampling(source.lab, clusters, IP.nClusters, IP.nSamples);

  otherwise
  disp('Invalid SAMPLE_METHOD');
end
samples.sourceSize = size(source.luminance);
toc;
if (OO.PLOT && ~IP.SUPERPIXEL)
  figure; imshow(source.image); title('Samples from source'); hold on;
  %Invert coordinates because it is a plot over an image.
  scatter(samples.idxs(2,:), samples.idxs(1,:), '.r'); hold off;

  figure(figs.ColorDist); title('Lab chrominance distribution (total x sampled)');
  scatter(samples.ab(1,:), samples.ab(2,:), 6, 'r');

  drawnow;
end

%% Feature extraction
dataName = IP.sourceFile(1:end-2-3);

try
  disp('Feature loading (try)'); tic;

  load(['./../temp/' dataName '_full']);
  toc;
catch
  disp('Feature extraction'); tic;

  [target_fv, target_fvl] = FeatureExtraction(target.luminance, FP);
  [samples_fv, samples_fvl] = FeatureExtraction(source.luminance, FP);
  toc;

  save(['./../temp/' dataName '_full'], 'target_fv', 'samples_fv', 'target_fvl', 'samples_fvl');
end

%Source Sampling
idxs = sub2ind(size(source.luminance), samples.idxs(1,:), samples.idxs(2,:));
%Structs receive values
samples.fv = samples_fv(:,idxs);
target.fv = target_fv;
samples.fvl = samples_fvl;
target.fvl = target_fvl;

%Clear structured variables
clear target_fv samples_fv target_fvl samples_fvl;

%% Superpixel extraction
if (IP.SUPERPIXEL)
  disp('Superpixel extraction'); tic;

  try
    disp('Superpixel loading'); tic
    
    load(['./../temp/' dataName '_sps']);
    toc;
  catch
    disp('Loading failed. Recomputing superpixels.'); tic
    
    [src_sp, src_lin_sp, src_centrs, src_nSP] = ...
      SuperpixelExtraction(source.luminance, IP.nSuperpixels, 'turbo');
    [tgt_sp, tgt_lin_sp, tgt_centrs, tgt_nSP] = ...
      SuperpixelExtraction(target.image, IP.nSuperpixels, 'turbo');
    toc;
    
    save(['./../temp/' dataName '_sps'], 'src_sp', 'tgt_sp', 'src_lin_sp', 'tgt_lin_sp', ...
      'src_centrs', 'tgt_centrs', 'src_nSP', 'tgt_nSP');
  end
  source.sp = src_sp; target.sp = tgt_sp;
  source.lin_sp = src_lin_sp; target.lin_sp = tgt_lin_sp;
  source.sp_centroids = src_centrs; target.sp_centroids = tgt_centrs;
  source.nSuperpixels = src_nSP; target.nSuperpixels = tgt_nSP;
  clear src_sp tgt_sp src_lin_sp tgt_lin_sp src_centrs tgt_centrs src_nSP tgt_nSP
  
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
end

%% Matching / Classification
disp('Feature matching / Classification in Feature Space'); tic;

SURF_mapping = true;

if (SURF_mapping)
  %% Source -> Target: initial label
  [~, nb_surf_idxs, nb_surf_dists ] = CombinedPDists(source.fv_sp, target.fv_sp, FP.feats1p);
  %Neighborhood classes from Source Automatic Labels
  nb_surf_classes = source.sp_clusters(nb_surf_idxs);
  
  %SURF initial classification 
  [surf_labels, ~, ~, ~, ~] = ...
    PredictSuperpixelsClassesKNN(nb_surf_classes, nb_surf_dists, IP.Kfs, IP.nClusters, clusters.mcCost);
  %Labels with the doubt flags.
%   surf_labels_final = surf_labels_c.*~surf_doubts + -1*surf_doubts;
  surf_labels_final = surf_labels;
  
  %Outputs:
  exp_labels{outLabels.p1, 1} = surf_labels_final;  
  exp_labels{outLabels.p1, 2} = 'Pass-1';
  if (OO.PLOT)
    surf_labels_img = CreateLabeledImage(surf_labels_final, target.sp, size(target.image));

    figure; imshow(surf_labels_img, []); colormap jet;
    title('Pass-1 labels');
    drawnow;
    clear surf_labels_img
  end
  
  clear nb_surf_classes nb_surf_dists nb_surf_idxs;
  
  %% Target -> Target: label refinement
  [~, nb_ld_idxs, nb_ld_dists] = CombinedPDists(target.fv_sp, target.fv_sp, FP.feats2p);
  nb_ld_classes = surf_labels_final(nb_ld_idxs);

  %Intra-image (remove the smallest distance which is the element itself).
  nb_ld_classes = [nb_ld_classes(:,2:end) nb_ld_classes(:,1)];
  nb_ld_idxs = [nb_ld_idxs(:,2:end) nb_ld_idxs(:,1)];
  nb_ld_dists = [nb_ld_dists(:,2:end) nb_ld_dists(:,1)];
%   nb_ld_dists(nb_ld_classes == -1) = Inf;
  
  %kNN majority
  exp_labels{outLabels.p2kNNm, 1} = modeTies(nb_ld_classes(:,1:IP.Kfs));
  exp_labels{outLabels.p2kNNm, 2} = 'Pass-2 kNNm';
  %kNN Equality
  [exp_labels{outLabels.p2kNNE, 1}, ...
   exp_labels{outLabels.p2kNNEcost, 1}, ...
   ~, final_scores, ~] = PredictSuperpixelsClassesKNN(nb_ld_classes, nb_ld_dists, IP.Kfs, IP.nClusters, clusters.mcCost);
  exp_labels{outLabels.p2kNNE, 2} = 'Pass-2 kNNE';
  exp_labels{outLabels.p2kNNEcost, 2} = 'Pass-2 cost kNNE';
  
  if (OO.PLOT)
    imClosestSP{outLabels.p2kNNm} = CreateLabeledImage(exp_labels{outLabels.p2kNNm,1}, target.sp, size(target.image));
    imClosestSP{outLabels.p2kNNE} = CreateLabeledImage(exp_labels{outLabels.p2kNNE,1}, target.sp, size(target.image));
    imClosestSP{outLabels.p2kNNEcost} = CreateLabeledImage(exp_labels{outLabels.p2kNNEcost,1}, target.sp, size(target.image));

    figure; imshow([imClosestSP{outLabels.p2kNNm}                 imClosestSP{outLabels.p2kNNE};
                    zeros(size(imClosestSP{outLabels.p2kNNm}))   imClosestSP{outLabels.p2kNNEcost}], []); colormap jet;
    title('Pass-2 labels: Scores> [kNNm - kNNE] ; Costs> [X - kNNE]');
    drawnow;
  end
  
  clear nb_ld_classes nb_ld_dists nb_ld_idxs;
end
  
%>Single-pass
[~, nb_s_idxs, nb_s_dists] = CombinedPDists(source.fv_sp, target.fv_sp, ...
  FP.feats1p + FP.feats2p);
nb_s_classes = source.sp_clusters(nb_s_idxs);
%kNN Equality
[~, exp_labels{outLabels.pSkNNE, 1}, ...
 ~, ~, ~] = PredictSuperpixelsClassesKNN(nb_s_classes, nb_s_dists, IP.Kfs, IP.nClusters, clusters.mcCost);
exp_labels{outLabels.pSkNNE, 2} = 'Single Pass';
%Matching
exp_labels{outLabels.pSmatch,1} = nb_s_idxs(:,1);
exp_labels{outLabels.pSmatch,2} = 'Single Matching idxs';

%>Color neighborhood
[~, neighbor_idxs, neighbor_dists] = CombinedPDists(source.fv_sp, target.fv_sp, FP.feats2p);
neighbor_classes = source.sp_clusters(neighbor_idxs);

clear nb_s_idxs nb_s_dists nb_s_classes

%% Edge-Aware Labeling/Relabeling
disp('Edge-Aware Labeling/Relabeling'); tic;
try
    load (['./../temp/' dataName '_eaclusters']);
    if(OO.PLOT)
      figure(73); imshow(eaClustersImg, []); colormap 'jet'
    end
catch
    [eaClusters, eaClustersImg] = EdgeAwareClustering(target);
    save (['./../temp/' dataName '_eaclusters'], 'eaClusters', 'eaClustersImg');
end

%Relabeling
exp_labels{outLabels.costkNNER,1} = EdgeAwareRelabeling2(eaClusters, exp_labels{outLabels.p2kNNEcost,1}, final_scores);
exp_labels{outLabels.costkNNER,2} = 'costkNNER';

% %Costs Labeling
% labelsEACPrE = EdgeAwareRelabeling(eaClusters, [], costsPrE);

if (OO.PLOT)
  imClosestSP{outLabels.costkNNER} = CreateLabeledImage(exp_labels{outLabels.costkNNER,1}, target.sp, size(target.image));

  
  figure; imshow(imClosestSP{outLabels.costkNNER}, []); colormap jet;
  title('Relabels> [costkNNE]');
  drawnow;
end

toc;

%% Output: Write the experiments images
disp('Experiment Outputs'); tic;

%Single x Double
for i = 1:4
  centroidImg = CreateCentroidImage(exp_labels{i,1},clusters.centroids,target.sp, target.image);
  imwrite(centroidImg, ['./../results/' batch_folder batch_out exp_labels{i,2} '.png'], 'png');
end

%Match:
matchImg = lab2rgb(CopyClosestSuperpixelAvgColor(source, target, exp_labels{outLabels.pSmatch,1}));
imwrite(matchImg, ['./../results/' batch_folder batch_out exp_labels{outLabels.pSmatch,2} '.png'], 'png');

%Classification
classImg = lab2rgb(CopyClosestSuperpixelFromClassAvgColor(source, target, ...
  neighbor_idxs, neighbor_classes, exp_labels{outLabels.p2kNNEcost,1}, 1));
imwrite(classImg, ['./../results/' batch_folder batch_out exp_labels{outLabels.p2kNNEcost,2} '.png'], 'png');

toc;

%% Color transfer:
disp('Color transfer + Save'); tic

%>Full method colorizations:
[tgt_scribbled, scribbles_mask] = CopyClosestSuperpixelFromClassAvgScribble(source, target, ...
      neighbor_idxs, neighbor_classes, exp_labels{outLabels.costkNNER,1}, IP.Kfs);
target.rgb = ColorPropagationLevin(lab2rgb(tgt_scribbled), target.luminance, scribbles_mask);
imwrite(target.rgb, ['./../results/' batch_folder batch_out dataName '_final_' exp_labels{outLabels.costkNNER,2} '.png'], 'png');

target.rgb_sat = IncreaseSaturation(target.rgb, 0.2);
imwrite(target.rgb_sat, ['./../results/' batch_folder batch_out dataName '_sat_' exp_labels{outLabels.costkNNER,2} '.png'], 'png');

toc;

%% Clear results 

clearvars -except batch_out batch_folder input_folder input_path list in_i;
