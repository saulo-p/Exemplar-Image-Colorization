function RGB = IncreaseSaturation(RGB, val)
  HSV = rgb2hsv(RGB);
  % "20% more" saturation:
  HSV(:, :, 2) = HSV(:, :, 2) * (1 + val);
  % or add:
  % HSV(:, :, 2) = HSV(:, :, 2) + 0.2;
  HSV(HSV > 1) = 1;  % Limit values
  RGB = hsv2rgb(HSV);
return;