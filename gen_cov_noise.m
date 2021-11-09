function noise_wav = gen_cov_noise(num_samps,num_chans,num_trials)


% 2 x 1 Gaussian noise vector with each variance. (orthogonal)
noise_wav = nan(num_samps,num_chans,num_trials);
for t=1:num_trials
    noiseVec = randn(num_chans,num_samps);
    c_cov = random('exp',3,num_chans,num_chans); % randomized covariance
    noiseVec = c_cov*noiseVec;
    
    noise_wav(:,:,t) = noiseVec';
end

figure(4); surf(cov(squeeze(noiseVec)')); shading interp; view(0,90); axis tight




