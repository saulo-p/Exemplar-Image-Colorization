%% -----------------------------------------------------------------
% Exemplar-based colorization algorithm
% Author: Saulo Pereira
%-------------------------------------------------------------------
clearvars -except batch_out batch_folder;
input_file = 'default';

%Structure to hold the outputs of the experiment
exp_labels = {};
exp_cols = {};
%Bloco 1.1
outLabels.kNNm = 1;
outLabels.kNNw = 2;
outLabels.kNNE = 3;
%Bloco 1.2
outLabels.costkNNw = 4;
outLabels.costkNNE = 5;
%Bloco 1.3
outColors.match = 1;
% outColors.costkNNw1 = 2;
% outColors.costkNNE1 = 3;
%Bloco 2
outLabels.costkNNwR = 6;
outLabels.costkNNER = 7;

%fieldsLabels = fieldnames(outLabels);
%Criar list@ das funcoes para executar em loop.

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

%% Luminance Remapping (source to target)
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
    toc;
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

%>Metric space construction
PDs = CombinedPDists(source.fv_sp, target.fv_sp, FP.featsWeights);
neighbor_idxs = zeros(size(PDs));
neighbor_dists = zeros(size(PDs));
for i = 1:size(PDs,1)
  [neighbor_dists(i,:), neighbor_idxs(i,:)] = sort(PDs(i,:));
end
clear PDs

%>Classification
neighbor_classes = source.sp_clusters(neighbor_idxs);
%kNN majority
exp_labels{outLabels.kNNm, 1} = modeTies(neighbor_classes(:,1:IP.Kfs));
exp_labels{outLabels.kNNm, 2} = 'kNNm';
%kNN weighted
[exp_labels{outLabels.kNNw, 1}, ...
 exp_labels{outLabels.costkNNw, 1}, ...
 ~, scoresW, costsW] = PredictSuperpixelsClassesKNN(neighbor_classes, neighbor_dists, -IP.Kfs, IP.nClusters, clusters.mcCost);
exp_labels{outLabels.kNNw, 2} = 'kNNw';
exp_labels{outLabels.costkNNw, 2} = 'costkNNw';
%kNN Equality
[exp_labels{outLabels.kNNE, 1}, ...
 exp_labels{outLabels.costkNNE, 1}, ...
 ~, scoresE, costsE] = PredictSuperpixelsClassesKNN(neighbor_classes, neighbor_dists, IP.Kfs, IP.nClusters, clusters.mcCost);
exp_labels{outLabels.kNNE, 2} = 'kNNE';
exp_labels{outLabels.costkNNE, 2} = 'costkNNE';

if (OO.PLOT)
  for i = 1:5
    imCentroids{i} = CreateCentroidImage(exp_labels{i,1}, clusters.centroids, target.sp, target.image*100);
  end

  figure; imshow([imCentroids{1}              imCentroids{2} imCentroids{3};
                  zeros(size(imCentroids{1})) imCentroids{4} imCentroids{5}], []); colormap jet;
  title('Locally assigned labels: Scores> [kNNm - kNNw - kNNE] ; Costs> [X - kNNw - kNNE]');
  drawnow;
end

toc;

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
exp_labels{outLabels.costkNNwR,1} = EdgeAwareRelabeling(eaClusters, exp_labels{outLabels.costkNNE,1}, []);
exp_labels{outLabels.costkNNwR,2} = 'costkNNwR';
exp_labels{outLabels.costkNNER,1} = EdgeAwareRelabeling(eaClusters, exp_labels{outLabels.costkNNE,1}, []);
exp_labels{outLabels.costkNNER,2} = 'costkNNER';

% %Costs Labeling
% labelsEACPrE = EdgeAwareRelabeling(eaClusters, [], costsPrE);

if (OO.PLOT)
  for i = 6:7
    imCentroids{i} = CreateCentroidImage(exp_labels{i,1}, clusters.centroids, target.sp, target.image*100);
  end
  
  figure; imshow([imCentroids{6} imCentroids{7} imCentroids{8}], []); colormap jet;
  title('Relabels> [kNNm - kNNw - kNNE]');
  
end

toc;

%% EXPERIMENTS: Write the experiments centroid images
for i = 1:7
  imCentroids{i} = CreateCentroidImage(exp_labels{i,1}, clusters.centroids, target.sp, target.image*100);
end

for i = 1:3
  imwrite(imCentroids{i}, ['./../results/' batch_folder batch_out dataName '_s11_' exp_labels{i,2} '_labels' '.png'], 'png');
end
for i = 4:5
  imwrite(imCentroids{i}, ['./../results/' batch_folder batch_out dataName '_s12_' exp_labels{i,2} '_labels' '.png'], 'png');
end

for i = 6:7
  imwrite(imCentroids{i}, ['./../results/' batch_folder batch_out dataName '_s2_' exp_labels{i,2} '_labels' '.png'], 'png');
end

%% Color transfer:
disp('Color transfer + Save'); tic

%>Matching colorization:
[exp_cols{outColors.match, 1}, ~] = knnsearch(source.fv_sp', target.fv_sp');
exp_cols{outColors.match,2} = 'match';
[tgt_scribbled, scribbles_mask] = CopyClosestSuperpixelAvgScribble(source, target, exp_cols{outColors.match, 1});
target.rgb = ColorPropagationLevin(lab2rgb(tgt_scribbled), target.luminance, scribbles_mask);
imwrite(target.rgb, ['./../results/' batch_folder batch_out dataName '_s13_' exp_cols{1,2} '.png'], 'png');

%>Classification colorizations: (using only the nearest neighbor color)
[tgt_scribbled, scribbles_mask] = CopyClosestSuperpixelFromClassScribble(source, target, ...
      neighbor_idxs, neighbor_classes, exp_labels{outLabels.costkNNw,1}, 1);
target.rgb = ColorPropagationLevin(lab2rgb(tgt_scribbled), target.luminance, scribbles_mask);
imwrite(target.rgb, ['./../results/' batch_folder batch_out dataName '_s13_' exp_labels{outLabels.costkNNw,2} '.png'], 'png');

[tgt_scribbled, scribbles_mask] = CopyClosestSuperpixelFromClassScribble(source, target, ...
      neighbor_idxs, neighbor_classes, exp_labels{outLabels.costkNNE,1}, 1);
target.rgb = ColorPropagationLevin(lab2rgb(tgt_scribbled), target.luminance, scribbles_mask);
imwrite(target.rgb, ['./../results/' batch_folder batch_out dataName '_s13_' exp_labels{outLabels.costkNNE,2} '.png'], 'png');

%>Full method colorizations:
%costKNNwR
[tgt_scribbled, scribbles_mask] = CopyClosestSuperpixelFromClassScribble(source, target, ...
      neighbor_idxs, neighbor_classes, exp_labels{outLabels.costkNNwR,1}, IP.Kfs);
target.rgb = ColorPropagationLevin(lab2rgb(tgt_scribbled), target.luminance, scribbles_mask);
imwrite(target.rgb, ['./../results/' batch_folder batch_out dataName '_final_' exp_labels{outLabels.costkNNwR,2} '.png'], 'png');
%costKNNER
[tgt_scribbled, scribbles_mask] = CopyClosestSuperpixelFromClassScribble(source, target, ...
      neighbor_idxs, neighbor_classes, exp_labels{outLabels.costkNNER,1}, IP.Kfs);
target.rgb = ColorPropagationLevin(lab2rgb(tgt_scribbled), target.luminance, scribbles_mask);
imwrite(target.rgb, ['./../results/' batch_folder batch_out dataName '_final_' exp_labels{outLabels.costkNNER,2} '.png'], 'png');

% for i = 2:length(exp_labels)
%   [tgt_scribbled, scribbles_mask] = CopyClosestSuperpixelFromClassScribble(source, target, ...
%       neighbor_idxs, neighbor_classes, exp_labels{i,1}, IP.Kfs);
%   tgt_scribbled = lab2rgb(tgt_scribbled);
%   target.rgb = ColorPropagationLevin(tgt_scribbled, target.luminance, scribbles_mask);
%   
%   %TEST:
% %   figure(1); imshow(target.rgb); title(img_gen{i,2});
% %   test_lab = rgb2lab(target.rgb);
% %   for c = 2:3
% %     figure(c); imshow(test_lab(:,:,c), []);
% %   end
%   
%   imwrite(target.rgb, ['./../results/' batch_folder batch_out dataName '_' exp_labels{i,2} '.png'], 'png');
% end

toc;