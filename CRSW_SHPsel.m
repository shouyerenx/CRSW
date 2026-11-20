function SHP = CRSW_SHPsel(slc, window, alpha, prefilter, workers)
% Select connected homogeneous pixels using complex covariance matrices.
[nlines, nwidths, nimages] = size(slc);
ntime = floor(nimages / 2);
if ntime < 2
    warning('CRSW:Selection', 'CCMmtf requires at least four acquisitions.');
end
if abs(alpha - 0.05) < 1e-4
    threshold = 9.4877;
elseif abs(alpha - 0.01) < 1e-4
    threshold = 13.2767;
else
    threshold = chi2inv(1 - alpha, 4);
end
C11 = zeros(nlines, nwidths, 'single');
C22 = zeros(nlines, nwidths, 'single');
C12 = zeros(nlines, nwidths, 'single') + 1j * zeros(nlines, nwidths, 'single');
for k = 1:ntime
    z1 = single(slc(:,:,2*k-1));
    z2 = single(slc(:,:,2*k));
    C11 = C11 + abs(z1).^2;
    C22 = C22 + abs(z2).^2;
    C12 = C12 + z1 .* conj(z2);
end
C11 = C11 + 1e-6;
C22 = C22 + 1e-6;
if prefilter
    filter_size = 3;
    kernel = fspecial('gaussian', [filter_size, filter_size], (filter_size - 1) / 6);
    C11 = imfilter(C11, kernel, 'replicate', 'conv');
    C22 = imfilter(C22, kernel, 'replicate', 'conv');
    C12 = imfilter(real(C12), kernel, 'replicate', 'conv') + ...
        1j * imfilter(imag(C12), kernel, 'replicate', 'conv');
    span = C11 + C22;
    average = fspecial('average', [filter_size, filter_size]);
    local_mean = imfilter(span, average, 'replicate');
    local_std = stdfilt(span, ones(filter_size, filter_size));
    variation = local_std ./ max(local_mean, eps('single'));
    enl = (0.52 ./ variation).^2;
    enl = max(enl, single(ntime));
    enl = min(enl, single(ntime * 25));
else
    enl = repmat(single(ntime), nlines, nwidths);
end
radius_row = floor(window(1) / 2);
radius_col = floor(window(2) / 2);
window_pixels = prod(window);
center = (2*radius_col+1)*radius_row + radius_col + 1;
C11_pad = padarray(C11, [radius_row, radius_col], NaN);
C22_pad = padarray(C22, [radius_row, radius_col], NaN);
C12_pad = padarray(C12, [radius_row, radius_col], NaN + 1j*NaN);
enl_pad = padarray(enl, [radius_row, radius_col], 'symmetric');
[pixel_index, line_index] = meshgrid(1:nwidths, 1:nlines);
index_matrix = sub2ind([nlines, nwidths], ...
    padarray(line_index, [radius_row, radius_col], 'symmetric'), ...
    padarray(pixel_index, [radius_row, radius_col], 'symmetric'));
total = nlines * nwidths;
pixel_sub = zeros(window_pixels, total, 'uint32');
brothers = zeros(total, 1, 'single');
log16 = log(single(16));
parfor (k = 1:total, workers)
    [row, col] = ind2sub([nlines, nwidths], k);
    row_pad = row + radius_row;
    col_pad = col + radius_col;
    row_range = row_pad-radius_row:row_pad+radius_row;
    col_range = col_pad-radius_col:col_pad+radius_col;
    a = reshape(C11_pad(row_range, col_range), window_pixels, 1);
    b = reshape(C22_pad(row_range, col_range), window_pixels, 1);
    c = reshape(C12_pad(row_range, col_range), window_pixels, 1);
    valid = ~isnan(a);
    if ~valid(center), continue; end
    a0 = a(center); b0 = b(center); c0 = c(center);
    enl_ref = enl_pad(row_pad, col_pad);
    rho = 1 - 7 / (4 * enl_ref);
    det_ref = max(a0*b0 - abs(c0)^2, eps('single'));
    det_window = max(a.*b - abs(c).^2, eps('single'));
    sum_a = a0+a; sum_b = b0+b; sum_c = c0+c;
    det_sum = max(sum_a.*sum_b - abs(sum_c).^2, eps('single'));
    statistic = -2 * rho * enl_ref * (log16 + log(det_ref) + log(det_window) - 2*log(det_sum));
    selected = false(window_pixels, 1);
    selected(valid) = statistic(valid) <= threshold;
    selected = reshape(selected, window(1), window(2));
    labels = bwlabel(selected, 8);
    center_label = labels(radius_row+1, radius_col+1);
    if center_label > 0
        selected = labels == center_label;
    else
        selected(:) = false;
        selected(center) = true;
    end
    selected = selected(:);
    indices = index_matrix(row_range, col_range);
    stored = zeros(window_pixels, 1, 'uint32');
    stored(selected) = indices(selected);
    pixel_sub(:,k) = stored;
    brothers(k) = sum(selected) - 1;
end
SHP.PixelSub = pixel_sub;
SHP.BroNum = reshape(brothers, [nlines, nwidths]);
SHP.CalWin = window;
end
