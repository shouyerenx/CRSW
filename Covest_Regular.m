function [coherence, sample_covariance] = Covest_Regular(slc, count)
% Estimate a regularized temporal coherence matrix.
[n, ~] = size(slc);
sample_covariance = (slc * slc') / count;
scale = sqrt(real(diag(sample_covariance)));
sample = sample_covariance ./ (scale .* scale' + 1e-12);
phase = angle(sample);
magnitude = abs(sample);

lags = zeros(n, 1);
lags(1) = 1;
diagonal_counts = 2 * (n - (0:n-1)');
diagonal_counts(1) = n;
for k = 1:n-1
    lags(k+1) = mean(diag(magnitude, k));
end
lags = max(lags, 0);
norm_sample = sum(magnitude(:).^2);
norm_target = sum(diagonal_counts .* (lags.^2));
distance = max(norm_sample - norm_target, 0);
if distance < 1e-10
    rho = 0;
else
    rho = ((1 / count) * norm_sample) / distance;
end
rho = min(max(rho, 0), 1);
regularized = (1 - rho) * magnitude + rho * toeplitz(lags);

log_bias = gammaln(count) + 0.5*log(pi) - log(2) - gammaln(count + 0.5);
bias = exp(log_bias);
noise_variance = max((1 / count) - bias^2, 1e-12);
smooth_lags = lags;
if n > 5
    smooth_lags(3:end) = movmean(lags(3:end), 5);
end
lag_index = (0:n-1)';
standard_error = sqrt(noise_variance ./ (n - lag_index));
margin = smooth_lags - bias;
scores = margin ./ (standard_error + 1e-9);
scores(1) = 100;
model_lags = max(margin, 1e-9);
model_lags(1) = 1;
target = toeplitz(model_lags);
distance = norm(regularized - target, 'fro')^2;
term = (1 - regularized.^2).^2;
term(1:n+1:end) = 0;
numerator = sum(term(:)) / (2 * count);
if distance < 1e-12
    global_weight = 0;
else
    global_weight = numerator / distance;
end
global_weight = max(0, min(1, global_weight));

offdiagonal = scores(2:end);
valid = isfinite(offdiagonal);
probability = zeros(size(offdiagonal));
if any(valid)
    low = min(offdiagonal(valid));
    high = max(offdiagonal(valid));
    if high - low > 0
        probability(valid) = (offdiagonal(valid) - low) / (high - low);
    end
end
probability = [1; max(0, min(1, probability))];
weights = global_weight .* probability + 1.0 * (1 - probability);
final_magnitude = (1 - toeplitz(weights)) .* regularized + toeplitz(weights) .* target;
combined = final_magnitude .* exp(1i * phase);
combined(1:n+1:end) = 1;
coherence = (combined + combined') / 2;
end
