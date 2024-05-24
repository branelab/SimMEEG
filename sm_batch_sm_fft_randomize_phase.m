function h = sm_batch_sm_fft_randomize_phase(h)
% use Alex Moiseev's program that will shuffle the phase using FFT of PCA components then put it back together as [sensors x samples x trials]

fprintf('Running FFT randomization of phases for sensor noise data.\n'); drawnow;


xdata = permute(h.sim_data.sens_noise(:,h.anatomy.sens.good_sensors,:),[2 1 3]);    % selecting only good channels
dims = size(xdata);
xdata = reshape(xdata,[dims(1) dims(2)*dims(3)]);

ydata = genSurrogateData(xdata);
h.sim_data.sens_noise(:,h.anatomy.sens.good_sensors,:) = permute(reshape(ydata,dims),[2 1 3]);
h.sim_data.sens_noise_scaled = h.sim_data.sens_noise./max(max(max(abs(h.sim_data.sens_noise)))); % scaled between -1 to 1
