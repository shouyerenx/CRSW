function [patch_ph, dates, shp_counts] = CRSW_run(cfg)
% Process one data set and save the phase cube and acquisition dates.
manifest = CRSW_manifest(cfg);
result_file = fullfile(char(cfg.output_dir), 'CRSW_result.mat');
assert(~isfile(result_file), 'CRSW:Output', 'Result already exists: %s', result_file);
if cfg.workers > 0
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('local', cfg.workers);
    elseif pool.NumWorkers < cfg.workers
        warning('CRSW:Workers', 'The existing pool has %d workers.', pool.NumWorkers);
    end
end
dates = manifest.dates;
nimages = numel(dates);
crop = manifest.crop;
height = crop(2)-crop(1)+1;
width = crop(4)-crop(3)+1;
stride = cfg.window_size - cfg.overlap;
nwindows = 1 + ceil(max(0, nimages-cfg.window_size) / stride);
patch_ph = zeros(height, width, nimages, 'single');
shp_counts = zeros(height, width, nwindows, 'single');
row_sections = section_limits(height, cfg.patches(1));
col_sections = section_limits(width, cfg.patches(2));
for row_patch = 1:size(row_sections,1)
    for col_patch = 1:size(col_sections,1)
        interior = [row_sections(row_patch,:), col_sections(col_patch,:)];
        halo = 2 * cfg.shp_half_window + 1;
        extended = [max(1,interior(1)-halo(1)), min(height,interior(2)+halo(1)), ...
                    max(1,interior(3)-halo(2)), min(width,interior(4)+halo(2))];
        bounds = extended + [crop(1)-1,crop(1)-1,crop(3)-1,crop(3)-1];
        fprintf('Processing patch %d/%d, %d/%d\n', row_patch, size(row_sections,1), ...
            col_patch, size(col_sections,1));
        [phase, counts] = process_patch(manifest, cfg, bounds);
        local_rows = (interior(1):interior(2)) - extended(1) + 1;
        local_cols = (interior(3):interior(4)) - extended(3) + 1;
        patch_ph(interior(1):interior(2),interior(3):interior(4),:) = ...
            single(phase(local_rows,local_cols,:));
        shp_counts(interior(1):interior(2),interior(3):interior(4),:) = ...
            single(counts(local_rows,local_cols,:));
    end
end
if ~isfolder(cfg.output_dir)
    [ok, message] = mkdir(cfg.output_dir);
    assert(ok, 'CRSW:Output', 'Cannot create output directory: %s', message);
end
assert(~isfile(result_file), 'CRSW:Output', 'Result already exists: %s', result_file);
temporary_file = [tempname(char(cfg.output_dir)), '.mat'];
cleanup = onCleanup(@() remove_temporary(temporary_file));
master_date = manifest.master_date;
save(temporary_file, 'patch_ph', 'dates', 'master_date', 'shp_counts', 'cfg', '-v7.3');
assert(~isfile(result_file), 'CRSW:Output', 'Result already exists: %s', result_file);
[ok, message] = movefile(temporary_file, result_file);
assert(ok, 'CRSW:Output', 'Cannot save result: %s', message);
fprintf('Saved %s\n', result_file);
end

function limits = section_limits(length_value, count)
section_size = ceil(length_value / count);
limits = zeros(count, 2);
for k = 1:count
    limits(k,1) = min(section_size*(k-1)+1, length_value-(count-k));
    limits(k,2) = min(section_size*k, length_value-(count-k));
end
end

function [phase, counts] = process_patch(manifest, cfg, bounds)
height = bounds(2)-bounds(1)+1;
width = bounds(4)-bounds(3)+1;
npix = height * width;
master_slc = CRSW_read_SARimage(manifest.slc{manifest.master_index}, cfg.nlines, ...
    manifest.width, cfg.slc_format, bounds);
master_diff = [];
if ~isempty(manifest.diff{manifest.master_index})
    master_diff = CRSW_read_SARimage(manifest.diff{manifest.master_index}, cfg.nlines, ...
        manifest.width, 'cpxfloat32', bounds);
end
total = numel(manifest.dates);
stride = cfg.window_size - cfg.overlap;
position = 0;
phase = [];
references = [];
counts = [];
anchor = [];
SHP = [];
while position < total
    indices = (position+1):min(position+cfg.window_size,total);
    if ~isempty(anchor) && numel(indices) <= cfg.overlap, break; end
    slc = CRSW_read_window(manifest, cfg, bounds, indices, master_slc, master_diff);
    if numel(indices) >= cfg.window_size || isempty(SHP)
        SHP = CRSW_SHPsel(slc, 2*cfg.shp_half_window+1, cfg.shp_alpha, ...
            logical(cfg.shp_prefilter), cfg.workers);
    end
    [current, counts, references, next_anchor, fused] = ...
        CRSW_core(slc, SHP, cfg, references, counts, anchor);
    if ~isempty(phase)
        n = size(fused, 2);
        tail = phase(:,end-n+1:end);
        valid = ~isnan(fused);
        tail(valid) = fused(valid);
        phase(:,end-n+1:end) = tail;
    end
    phase = [phase, current];
    anchor = next_anchor;
    position = position + stride;
end
assert(size(phase,2) == total, 'CRSW:Processing', 'Incomplete phase history.');
if size(references,2) > 1
    reference_phase = CRSW_reference(references, SHP, cfg, [height,width]);
    expanded = zeros(size(phase), 'single');
    current_col = 0;
    for k = 1:size(references,2)
        n = stride;
        if k == size(references,2), n = size(phase,2)-current_col; end
        if n <= 0, break; end
        expanded(:,current_col+1:current_col+n) = repmat(reference_phase(:,k), 1, n);
        current_col = current_col+n;
    end
    phase = angle(exp(1j * phase) .* exp(1j * expanded));
end
phase = reshape(phase, height, width, total);
counts = reshape(counts, height, width, size(counts,2));
end

function remove_temporary(filename)
if isfile(filename), delete(filename); end
end
