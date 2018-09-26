function [im_centroid] = CreateCentroidImage(sp_labels, cl_centroids, image_sps, img_gray)
%

img_size = size(img_gray);
im_centroid = zeros(img_size(1), img_size(2), 3);
im_centroid(:,:,1) = img_gray;

for i = 1:length(sp_labels)
  if (sp_labels(i) == -1)
    continue;
  end
  
  mask = (image_sps == i);
  for c = 2:3  
    im_centroid(:,:,c) = im_centroid(:,:,c) + cl_centroids(sp_labels(i),c-1)*mask;
  end
end

im_centroid = lab2rgb(im_centroid);
im_centroid(1,1) = -1; %For comparison with classification

end

