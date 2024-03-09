%% BANDWIDTH PLOT
clear all

filenames=[dir('*43838*.h5')];
Fignum=1096;

remove_on=1;
remove_threshold=40;

num_ca=1;
signal_per_ion=90;
noise=22.5;


%% PREPARE DATA
% set count thresholds
threshold=signal_per_ion/log(1+signal_per_ion/noise);
threshold2=180;sqrt(2)*signal_per_ion;
if num_ca==1
    threshold2=sqrt(2)*signal_per_ion*10;
end

data=[];
% 1st: AOM freq; 2nd: LIF; ; 3th: carrier freq; 4th: secular frequency

for i=1:max(size(filenames))
    clear data1
    data1 = h5read(filenames(i).name,'/datasets/results');
    data=[data; data1'];

    % get experimental parameters
    eggs_ht=h5readatt(filenames(i).name,'/arguments/','time_eggs_heating_ms');  %ms
    eggs_att=h5readatt(filenames(i).name,'/arguments/','att_eggs_heating_db');  %dB
    time=h5readatt(filenames(i).name,'/arguments/','time_sideband_readout_us');  %dB

    eggs_freq_carrier_mhz = h5readatt(filenames(i).name,'/arguments/','freq_eggs_heating_carrier_mhz_list');
    eggs_ampl_carrier_pct = h5readatt(filenames(i).name,'/arguments/','ampl_eggs_dynamical_decoupling_pct');
    eggs_ampl_rsb_pct = h5readatt(filenames(i).name,'/arguments/','ampl_eggs_heating_rsb_pct');
    eggs_ampl_bsb_pct = h5readatt(filenames(i).name,'/arguments/','ampl_eggs_heating_bsb_pct');
    phase_rsb_turns = h5readatt(filenames(i).name,'/arguments/','phase_eggs_heating_rsb_turns_list');
    phase_bsb_turns = h5readatt(filenames(i).name,'/arguments/','phase_eggs_heating_bsb_turns');

end

figure;plot(data(:,2))
if remove_on==1
    data=remove_zero(data,remove_threshold);
end


%% COLLATE DATA
AOM_freq=sort(unique(data(:,1)));


%% PROCESS DATA
%threshold counts
counts_thresholded = ((data(:, 2) > threshold) + (data(:, 2) > threshold2)) / num_ca;

% separate rsb and bsb
rsb_data = counts_thresholded(data(:,1)==AOM_freq(1));
bsb_data = counts_thresholded(data(:,1)==AOM_freq(2));

% get sideband rabi ampl.
c = 1. - mean(rsb_data);
rsb_ampl = 1. - mean(rsb_data);
bsb_ampl = 1. - mean(bsb_data);
err_rsb = std(rsb_data)/sqrt(length(rsb_data));
err_bsb = std(bsb_data)/sqrt(length(rsb_data));
sprintf('RSB: %.3f +/- %.3f \t BSB: %.3f +/- %.3f', rsb_ampl, err_rsb, bsb_ampl, err_bsb)

% calculate phonon number
ratio_sb = rsb_ampl / bsb_ampl;
err_ratio_sb = ratio_sb * sqrt((err_rsb/rsb_ampl)^2. + (err_bsb/bsb_ampl)^2.);
n_coherent = coherent_c(ratio_sb);
err_n_coherent = (coherent_c(ratio_sb + 2 * err_ratio_sb) - coherent_c(ratio_sb - 2 * err_ratio_sb))/4.;
n_ratio = ratio_sb / (1. - ratio_sb);
err_n_ratio = n_ratio * bsb_ampl * sqrt((err_rsb / rsb_ampl)^2. + (err_bsb / bsb_ampl)^2.);
n_coherent = coherent_c(n_ratio);

% PRINT RESULT
sprintf('<n> = %.3f +/- %.3f', n_coherent, err_n_coherent)
% sprintf('<n> = %.3f +/- %.3f\nRSB Phase: %.3f\nBSB Phase: %.3f', n_coherent, err_n_coherent, phase_rsb_turns, phase_bsb_turns)
