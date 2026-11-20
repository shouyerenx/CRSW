load('./CRSW_output/CRSW_result.mat')
m =45;
n =85;
figure
phase_m_my = patch_ph(:, :, m); 
phase_n_my = patch_ph(:, :, n); 
Phase_my_recon = angle(exp(1j * (phase_n_my - phase_m_my))); 
imagesc(Phase_my_recon, [-pi, pi]);
axis tight;
set(gca, 'XTick', [], 'YTick', []);
axis image;
colormap(jet);