function [f_m,f_med,f_ptil_h,f_ptil_l] = stats_ptc(f,prc_h,prc_l)

% The first two dims of f must be space

f_m = squeeze(nanmean(f,[1 2]));
f_med = squeeze(prctile(f,50,[1 2]));
f_ptil_h = squeeze(prctile(f,prc_h,[1 2]));
f_ptil_l = squeeze(prctile(f,prc_l,[1 2]));

end
