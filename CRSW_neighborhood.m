function [valid, mask, center, row_index, col_index] = CRSW_neighborhood(SHP, pixel, half_window, image_size, min_shp)
% Return the selected neighborhood and its wrapped raster indices.
height = image_size(1); width = image_size(2);
row = mod(pixel-1, height) + 1;
col = ceil(pixel / height);
mask = logical(reshape(SHP.PixelSub(:,pixel), 2*half_window(1)+1, 2*half_window(2)+1));
valid = sum(mask(:)) >= min_shp;
row_index = rem(height + rem((row-half_window(1):row+half_window(1))-1, height), height) + 1;
col_index = rem(width + rem((col-half_window(2):col+half_window(2))-1, width), width) + 1;
center = sub2ind(size(mask), half_window(1)+1, half_window(2)+1);
end
