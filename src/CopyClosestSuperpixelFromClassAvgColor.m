function lab_out = CopyClosestSuperpixelFromClassAvgColor(source, target, ...
  neighbor_idxs, neighbor_classes, labels, kCT)
%Each superpixel receives the color of the closest superpixel from the
%majority class.

%Output image
lab_out = zeros([size(target.image) 3]);
lab_out(:,:,1) = target.luminance*100;

for i = 1:target.nSuperpixels
  if (labels(i) == -1)
    tgt_mask = (target.sp==i);
    lab_out(:,:,1) = lab_out(:,:,1).*~tgt_mask;
    %if label is marked as doubt, assign black (for debug purposes)  
    continue
  end

  % Instances from chosen class
  [~, majority_instances] = find(neighbor_classes(i,:) == labels(i));

  %Matching superpixels ROI masks
  tgt_mask = (target.sp==i);
  src_mask = zeros(size(source.sp));
  for ni = 1:min([kCT length(majority_instances)])
    src_mask = src_mask | (source.sp==neighbor_idxs(i,majority_instances(ni)));
  end
  
  for c = 2:3
    %Prototype color transfer (Superpixel average)
    mask_c = source.lab(:,:,c).*src_mask;
    avg_sp = sum(mask_c(:))/length(find(src_mask));
    lab_out(:,:,c) = lab_out(:,:,c) + avg_sp*tgt_mask;
  end
end

end

