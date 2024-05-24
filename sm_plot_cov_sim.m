function sm_plot_cov_sim(varargin)
global h


figure(200+h.menu_synthetic_noise_cov_type.Value); clf; 
ax = subplot_axes(3,2,.1,.075, 0, 0, 0); 


%% calculating covariances
dims = size(h.sim_data.sens_noise); 
noise_r = nan(dims(2),dims(2),dims(3));  signal_r = noise_r; final_r = noise_r;
for t=1:size(h.sim_data.sens_noise,3)
    noise_r(:,:,t) = cov(squeeze(h.sim_data.sens_noise_final(h.cfg.study.base_samps,:,t)));
    signal_r(:,:,t) = cov(squeeze(h.sim_data.sens_sig_data(h.cfg.study.act_samps,:,t)));
    final_r(:,:,t) = cov(squeeze(h.sim_data.sens_final(h.cfg.study.act_samps,:,t)));
end
noise_r = squeeze(nanmean(noise_r,3));  % averaged trial-by-trial covariances
signal_r = squeeze(nanmean(signal_r,3)); 
final_r = squeeze(nanmean(final_r,3)); 

noise_revk = cov(squeeze(nanmean(h.sim_data.sens_noise(h.cfg.study.base_samps,:,:),3)));    % evoked covariance
signal_revk = cov(squeeze(nanmean(h.sim_data.sens_sig_data(h.cfg.study.act_samps,:,:),3)));
final_revk = cov(squeeze(nanmean(h.sim_data.sens_final(h.cfg.study.act_samps,:,:),3)));

%% getting starting scales
noise_r_scale = quantile(abs(reshape(noise_r,numel(noise_r),1)),.90);  
noise_revk_scale = quantile(abs(reshape(noise_revk,numel(noise_revk),1)),.90);  
signal_r_scale = quantile(abs(reshape(signal_r,numel(signal_r),1)),.90);  
signal_revk_scale = quantile(abs(reshape(signal_revk,numel(signal_revk),1)),.90);  
final_r_scale = quantile(abs(reshape(final_r,numel(final_r),1)),.90);  
final_revk_scale = quantile(abs(reshape(final_revk,numel(final_revk),1)),.90);  

%% Plotting Covariances 
% Noise averaged trial
axes(ax(1)); cla;
surf(squeeze(nanmean(noise_r,3))); view(0,90); shading interp; axis tight; caxis([-noise_r_scale noise_r_scale]);
title('Noise Covariance (averaged-trial)'); colorbar; colormap(jet);
% Noise evoked
axes(ax(2)); cla;
surf(squeeze(nanmean(noise_revk,3))); view(0,90); shading interp; axis tight; caxis([-noise_revk_scale noise_revk_scale]);
title('Noise Covariance (evoked)'); colorbar; colormap(jet);
% Signal averaged trial
axes(ax(3)); cla;
surf(squeeze(nanmean(signal_r,3))); view(0,90); shading interp; axis tight; caxis([-signal_r_scale signal_r_scale]);
title('Signal Covariance (averaged-trial)'); colorbar; colormap(jet);
% Signal evoked
axes(ax(4)); cla;
surf(squeeze(nanmean(signal_revk,3))); view(0,90); shading interp; axis tight; caxis([-signal_revk_scale signal_revk_scale]);
title('Signal Covariance (evoked)'); colorbar; colormap(jet);
% Final averaged trial
axes(ax(5)); cla;
surf(squeeze(nanmean(final_r,3))); view(0,90); shading interp; axis tight; caxis([-final_r_scale final_r_scale]);
title('Final Covariance (averaged-trial)'); colorbar; xlabel('Channel #'); ylabel('Channel #');
% Final evoked
axes(ax(6)); cla;
surf(squeeze(nanmean(final_revk,3))); view(0,90); shading interp; axis tight; caxis([-final_revk_scale final_revk_scale]);
title('Final Covariance (evoked)'); colorbar; colormap(jet);




