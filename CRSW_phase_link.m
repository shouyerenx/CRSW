function [phasor, phase] = CRSW_phase_link(coherence)
% Estimate the phase history from a temporal coherence matrix.
gamma = abs(coherence);
[~, eigenvalues] = eig(double(gamma));
count = 1;
while min(diag(eigenvalues)) < 0
    diagonal = diag(gamma);
    gamma = 1e-6 * gamma;
    gamma(logical(eye(size(gamma)))) = diagonal;
    [~, eigenvalues] = eig(double(gamma));
    count = count + 1;
    if count > size(coherence, 1), break; end
end
[vectors, ~] = svd(pinv(gamma) .* coherence);
phasor = vectors(:, end);
phasor = sqrt(size(coherence, 1)) * phasor / norm(phasor);
phase = angle(phasor);
phasor = exp(1j * (phase - phase(1)));
phase = angle(phasor);
end
