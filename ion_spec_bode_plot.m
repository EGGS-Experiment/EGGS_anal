%% ISA: PHASE SWEEP & WSEC SWEEP W/SINGLE POINT SIDEBAND READOUT
%Use this code when you Ion Spectrum Analyzer
%with single sideband values, a single readout time, multiple secular frequency values and one
%detuning i.e. on resonance.

% clear
%figure;plot(data(:,2))
%figure;histogram(data(:,2))
data=[];

%% CONFIGURE
date_path = '\\eric.physics.ucla.edu\groups\motion\Data\2024-03\2024-03-08'; 
rid_dj = '51120';

rid_str = strcat('*', rid_dj, '*.h5');
filenames=[dir(fullfile(date_path, rid_str))];
Fignum=1098;

% REMOVE ZERO
remove_on=0;
remove_threshold=100;

% FLUORESCENCE
num_ca=1;
noise=25;
signal_per_ion=160;

% CALCULATE THRESHOLDS
threshold=signal_per_ion/log(1+signal_per_ion/noise);
threshold2=180;sqrt(2)*signal_per_ion;
if num_ca==1
    threshold2=sqrt(2)*signal_per_ion*10;
end


%% IMPORT DATA

for i=1:max(size(filenames))
    clear data1
    file = fullfile(date_path, filenames(i).name); 
    data1 = h5read(file,'/datasets/results');
    data=[data; data1'];

    % get experimental parameters
    eggs_ht=h5readatt(file,'/arguments/','time_eggs_heating_ms');  %ms
    %eggs_ht = 1.0;
    eggs_att=h5readatt(file,'/arguments/','att_eggs_heating_db');  %dB
    time=h5readatt(file,'/arguments/','time_sideband_readout_us');  %dB

    eggs_freq_carrier_mhz = h5readatt(file,'/arguments/','freq_eggs_heating_carrier_mhz_list');
    eggs_ampl_carrier_pct = h5readatt(file,'/arguments/','ampl_eggs_dynamical_decoupling_pct');
    eggs_ampl_rsb_pct = h5readatt(file,'/arguments/','ampl_eggs_heating_rsb_pct');
    eggs_ampl_bsb_pct = h5readatt(file,'/arguments/','ampl_eggs_heating_bsb_pct');
end

% figure;plot(data(:,2))
if remove_on==1
    data=remove_zero(data,remove_threshold);
end


%% PRE-PROCESS DATASET
% separate data by column
% 1st: AOM freq; 2nd: LIF; ; 3th: carrier freq; 4th: secular frequency;

aom_freq_index =        1;
carrier_freq_index =    3;
secular_freq_index =    4;
rsb_phase_index =       5;
AOM_freq=unique(data(:,aom_freq_index));
carrier_freq_all=data(:, carrier_freq_index);
secular_freq_all=data(:, secular_freq_index);
rsb_phase_all=data(:, rsb_phase_index);

% get sweep values and range
carrier_freq=unique(carrier_freq_all);
secular_freq=unique(secular_freq_all);
rsb_phase=unique(rsb_phase_all);

num_AOM=max(size(AOM_freq)) ;
num_carrier=max(size(carrier_freq));
num_secular=max(size(secular_freq));
num_bsb=max(size(rsb_phase));

phase=[];
all_LIF_b=[];
all_LIF_r=[];
all_fits_phonon=[];
err_phonon=[];

% guess parameters for separating data
freq2 = 2*(AOM_freq)*2.32830644e-7;
guess_carrier_freq_mhz = mean(freq2);
guess_secular_res = mean(secular_freq) / 1e6;
%guess_secular_res = 0.7579;


%% PLOT AND FIT DATA
% in this code, secular and carrier are flipped to quickly plot data
for j=1:num_bsb
    for m=1: num_secular
        clear LIF
        LIF=[];
        for i=1: num_AOM
                index=data(:, aom_freq_index)==AOM_freq(i) ...
                    &data(:, secular_freq_index)==secular_freq(m)...
                    &data(:, rsb_phase_index)==rsb_phase(j);

               order1=data(index,2)>threshold &data(index ,2)<threshold2;
               order2=data(index ,2)>threshold2;
               LIF(i)=mean(order1+2*order2)/num_ca;   
               LIF_err(i)=std(order1+2*order2)/num_ca/sqrt(max(size(order1)));
        end
           
        freq_r=freq2(freq2<guess_carrier_freq_mhz);
        LIF_r=(LIF(freq2<guess_carrier_freq_mhz));
        LIF_r_err= LIF_err(freq2<guess_carrier_freq_mhz);

        freq_b=freq2(freq2>guess_carrier_freq_mhz);
        LIF_b=(LIF(freq2>guess_carrier_freq_mhz));
        LIF_b_err= LIF_err(freq2>guess_carrier_freq_mhz);
             
        all_LIF_r(m)=1-LIF_r;
        all_LIF_b(m)=1-LIF_b;
    end
    
    % PLOT RSB
    figure(j+2*Fignum);subplot(2,2,1);hold on
    plot(secular_freq*1e-6,all_LIF_r','or')
    xlabel('Secular Freq. (MHz)');ylabel('D state');ylim([0, 1]);
    title(sprintf('%.2f', rsb_phase(j)));
    % PLOT BSB
    figure(j+2*Fignum);subplot(2,2,2);hold on
    plot(secular_freq*1e-6,all_LIF_b','ob')
    xlabel('Secular Freq. (MHz)');ylabel('D state');ylim([0, 1]);

    % RSB/BSB RATIO TO DISPLACEMENT
    % add 0.001 to prevent NaN error from forming
    all_LIF_b = all_LIF_b + 0.0001;
    all_ratio=all_LIF_r./all_LIF_b;    
    all_ratio(all_ratio<0)=0;
    all_ratio(all_ratio>0.791976610483122)=0.791976610483122;
    [all_phonon,all_err]=coherent_c(all_ratio);
    %all_phonon=squeezing_s(all_ratio);
    % PLOT RSB/BSB RATIO
    figure(j+2*Fignum);subplot(2,2,3);hold on
    plot(secular_freq*1e-6,all_ratio,'o')
    title('RSB/BSB Ratio')
    center=mean( secular_freq(all_ratio==max(all_ratio))*1e-6 );

    % fit <n> with sinc^2
    % guess parameters
    guess_phonon_offset = mean([min(all_phonon'), median(all_phonon')]);
    guess_phonon_ampl = (1 - min(all_phonon')) - guess_phonon_offset;
    guess_phonon_ampl = 2 * asin(sqrt(guess_phonon_ampl)) / (pi * time);
    fo_phonon = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, guess_secular_res-0.005, 1e3*eggs_ht, 0],...
               'Upper',[2*pi*20e-4, guess_secular_res+0.005, 1e3*eggs_ht, 0.200],...
               'StartPoint',[guess_phonon_ampl, guess_secular_res, 1e3*eggs_ht, guess_phonon_offset]);
    ftphonon = fittype('2*a^2/(x-b+1e-10)^2*(sin(2*pi*(x-b+1e-10)*c/2))^2+d','options',fo_phonon);
    [curvephonon,gofphonon] = fit(secular_freq*1e-6,all_phonon',ftphonon);
    all_fits_phonon(j) = curvephonon.b;

    % store fit results
    phonon_all(j)=curvephonon(curvephonon.b)-curvephonon.d;
    %phonon_all(j)=mean(all_phonon);
    phase(j)=rsb_phase(j);

    % CALCULATE ERROR
    all_std=confint(curvephonon);
    err_a=(all_std(2,1)-all_std(1,1))/4 ;

    if isnan((all_std(2,3)-all_std(1,3))/4)==1
        err_c=0;
    else
        err_c=(all_std(2,3)-all_std(1,3))/4;
        err_c=0.04;
    end

    if isnan((all_std(2,4)-all_std(1,4))/4)==1
        err_d=0;
    else
        err_d=(all_std(2,4)-all_std(1,4))/4;
    end
    
    % PLOT

    last_file = fullfile(date_path, filenames(end).name); 
    figure(j+2*Fignum);title([last_file(67:71)]);
    subplot(2,2,4);hold on
    plot(secular_freq*1e-6,all_phonon,'o');plot(curvephonon)
    legend('hide');title('phonon');xlabel('Secular Freq. (MHz)');ylabel('<n>') 
end


%% SUMMARY PLOT
% fit data
fo_phase= fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, 0],...
               'Upper',[0.9, 2.*pi],...
               'StartPoint',[0.25, 1.*pi]);
ftphase = fittype('a*sin(2*pi*x+c)+a','options',fo_phase);
[curvephase,gofphase] = fit(phase',phonon_all',ftphase);

% plot!
figure(10000);plot(phase,phonon_all,'o')
hold on;plot(curvephase);legend('hide');
xlabel('Phase (turns)');ylabel('<n>');
% xlabel('ISA Carrier Freq. (10x MHz)');

% calculate error
err_phasefit = confint(curvephase);
err_ampl = 2. * (err_phasefit(2,1) - err_phasefit(1,1)) / 4.;
err_phase = (err_phasefit(2,2) - err_phasefit(1,2)) / 4.;

% print fit values
params_fit1 = sprintf('Ampl:  %.4f +/- %.4f', curvephase.a*2., err_ampl);
params_fit2 = sprintf('Phase:  %.4f +/- %.4f', 2.*.5*(0.25-curvephase.c/(2.*pi)), err_phase/(2.*pi));
% params_fit2 = sprintf('Phase:  %.4f +/- %.4f', (curvephase.c/(2.*pi)), err_phase/(2.*pi));
text(0.3, 0.1, params_fit1)
text(0.3, 0.15, params_fit2)

% print experimental parameters

last_file = fullfile(date_path, filenames(end).name); 
params_exp = sprintf('RID: %s    %.1f MHz\n %.1f dB: %.1f/%.1f/%.1f', ...
    last_file(67:71), eggs_freq_carrier_mhz, eggs_att, eggs_ampl_rsb_pct, eggs_ampl_bsb_pct, eggs_ampl_carrier_pct);
% params_exp = sprintf('RID: %s  \n %.1f dB: %.1f/%.1f/%.1f', ...
%     last_file(67:71), eggs_att, eggs_ampl_rsb_pct, eggs_ampl_bsb_pct, eggs_ampl_carrier_pct);
title(params_exp);
sprintf('Opt. Carr. Ampl: %.2f%%\nwsec: %.3f', sqrt(0.88/(curvephase.a*2.)) * eggs_ampl_carrier_pct, mean(all_fits_phonon)*1000)
