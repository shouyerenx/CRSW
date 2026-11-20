function correction = CRSW_reference(references, SHP, cfg, image_size)
% Link all compressed references for the current spatial patch.
npix = prod(image_size);
nref = size(references, 2);
stack = reshape(references, image_size(1), image_size(2), nref);
correction = zeros(npix, nref, 'single');
half = cfg.shp_half_window;
parfor (pixel = 1:npix, cfg.workers)
    [valid, mask, ~, rr, cc] = CRSW_neighborhood(SHP, pixel, half, image_size, cfg.min_shp);
    if ~valid, continue; end
    block = stack(rr(:), cc(:), :);
    block = reshape(block, prod(2*half+1), nref);
    block = block(mask(:), :);
    block = block(~any(isnan(block), 2), :).';
    count = size(block, 2);
    if count < 1, continue; end
    coherence = Covest_Regular(block, count);
    [~, local] = CRSW_phase_link(coherence);
    correction(pixel,:) = local.';
end
end
