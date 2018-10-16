function remapped = luminance_remap(source, target, AUTO_COL)
    % assumes that the source is Lab with luminance channel (0-100)
    % assumes that target is a grayscale image (0-1)
    
    if(AUTO_COL)
        disp('Auto-colorization');
        remapped = target;
        return;
    end
    
    source_luminance = source(:,:,1)/100;
    
    mu_a = mean(source_luminance(:));
    mu_b = mean(target(:));
    sigma_a = std(source_luminance(:));
    sigma_b = std(target(:));
    
    remapped = (sigma_b/sigma_a)*(source_luminance - mu_a) + mu_b;
    %% Test code;
%     matched = imhistmatch(source_luminance, target);
%     
%     figure;
%     subplot(4,1,1); histogram(source_luminance); title('Original histogram');
%     subplot(4,1,2); histogram(target(:)); title('Target histogram');
%     subplot(4,1,3); histogram(remapped(:)); title('Linearly remapped histogram');
%     subplot(4,1,4); histogram(matched(:)); title('Matched histogram');
%     
%     figure;
%     imshow([remapped matched]);
%     title('Linear remap (left) and Match (right)');
%     
%     remapped = matched;
end