function [phase, counts, references, next_anchor, fused] = CRSW_core(slc, SHP, cfg, references, counts, anchor)
% Link one temporal window and update its compressed reference.
[height, width, nimages] = size(slc);
npix = height * width;
first = isempty(references);
if first, start = 1; else, start = cfg.overlap + 1; end
nout = max(0, nimages-start+1);
phase = zeros(npix, nout);
next_anchor = nan(npix, cfg.overlap, 'single');
fused = nan(npix, cfg.overlap, 'single');
new_reference = zeros(npix, 1, 'like', 1i);
new_counts = zeros(npix, 1);
half = cfg.shp_half_window;
parfor (pixel = 1:npix, cfg.workers)
    [valid, mask, center, rr, cc] = CRSW_neighborhood(SHP, pixel, half, [height,width], cfg.min_shp);
    if ~valid, continue; end
    block = slc(rr(:), cc(:), :);
    block = reshape(block, prod(2*half+1), nimages);
    keep = ~any(isnan(block), 2);
    block = block(keep, :).';
    if isempty(block) || ~keep(center), continue; end
    center_real = sum(keep(1:center));
    selected = mask(:);
    selected = selected(keep);
    selected_block = block(:,selected);
    count = size(selected_block, 2);
    if count < 1, continue; end
    coherence = Covest_Regular(selected_block, count);
    [phasor, local] = CRSW_phase_link(coherence);
    delta = 0;
    row_fused = nan(1, cfg.overlap, 'single');
    if ~first && ~isempty(anchor)
        n_overlap = min([cfg.overlap, size(anchor,2), numel(local)]);
        if n_overlap > 0
            previous = anchor(pixel,end-n_overlap+1:end);
            if ~any(isnan(previous))
                current = local(1:n_overlap).';
                difference = angle(exp(1j * (previous - current)));
                alignment = angle(mean(exp(1j * difference)));
                for iteration = 1:2
                    residual = angle(exp(1j * (difference - alignment)));
                    weights = 1 ./ (1 + (residual / 0.8).^2);
                    alignment = angle(sum(weights .* exp(1j * difference)));
                end
                aligned = mod(current + alignment + pi, 2*pi) - pi;
                fused_value = angle(exp(1j * previous) + exp(1j * aligned));
                refined = angle(exp(1j * (fused_value - current)));
                delta = angle(mean(exp(1j * refined)));
                row_fused(end-n_overlap+1:end) = fused_value;
            end
        end
    end
    fused(pixel,:) = row_fused;
    rotated = mod(local + delta + pi, 2*pi) - pi;
    if nout > 0, phase(pixel,:) = rotated(start:end).'; end
    row_anchor = nan(1, cfg.overlap, 'single');
    nsave = min(cfg.overlap, numel(rotated));
    row_anchor(end-nsave+1:end) = rotated(end-nsave+1:end).';
    next_anchor(pixel,:) = row_anchor;
    phasor = phasor * exp(1j * delta);
    phasor = phasor / max(1, norm(phasor));
    compressed = phasor' * block;
    new_reference(pixel) = compressed(center_real);
    new_counts(pixel) = count;
end
references = [references, new_reference];
counts = [counts, new_counts];
end
