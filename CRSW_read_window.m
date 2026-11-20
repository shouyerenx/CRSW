function slc = CRSW_read_window(manifest, cfg, bounds, indices, master_slc, master_diff)
% Remove the reference topographic phase from each SLC.
slc = complex(zeros(bounds(2)-bounds(1)+1, bounds(4)-bounds(3)+1, numel(indices), 'single'));
for k = 1:numel(indices)
    index = indices(k);
    raw = CRSW_read_SARimage(manifest.slc{index}, cfg.nlines, manifest.width, cfg.slc_format, bounds);
    if index == manifest.master_index
        if isempty(master_diff)
            corrected_diff = ones(size(raw), 'single');
        else
            corrected_diff = master_diff .* 0 + (1 + 0i);
        end
    else
        corrected_diff = CRSW_read_SARimage(manifest.diff{index}, cfg.nlines, manifest.width, ...
            'cpxfloat32', bounds);
        corrected_diff = corrected_diff ./ abs(corrected_diff);
    end
    original_diff = exp(1i * angle(master_slc .* conj(raw)));
    slc(:,:,k) = raw .* (original_diff .* conj(corrected_diff));
end
end
