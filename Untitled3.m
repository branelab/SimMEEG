



peak_tfr = h.current_inv_peak_fc_data;
true_tfr = h.current_inv_true_fc_data;
act_samps = h.inv_soln(h.current_inv_soln).params.act_samps;
peak_tfr(abs(peak_tfr)>0)=10;   % values of 10 = significant index; 0 = not significant 
true_tfr(abs(true_tfr)>0)=5;   % values of 5 = significant index; 0 = not significant 


diff_tfr = peak_tfr-true_tfr; 
hit_tfr = nan(size(diff_tfr)); hit_tfr(diff_tfr==5)=1;
miss_tfr = nan(size(diff_tfr)); miss_tfr(diff_tfr== -5) = 1;
fa_tfr = nan(size(diff_tfr)); fa_tfr(diff_tfr== 10) = 1;
cr_tfr = nan(size(diff_tfr)); cr_tfr(diff_tfr== 0) = 1;


num_peak_sig = squeeze(nansum(nansum(peak_tfr(:,:,act_samps),3),1));
num_true_sig = squeeze(nansum(nansum(true_tfr(:,:,act_samps),3),1));


