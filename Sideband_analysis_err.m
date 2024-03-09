%% SIDEBAND_ANALYSIS.M: PROCESS A SIDEBAND SPECTRUM READOUT TO MEASURE <n>
%figure;plot(data(:,2))
%figure;histogram(data(:,2))
clear;

% file/display
filenames=dir(['*44513*.h5']);
Fignum=271;

% remove data consisting of consecutive 0s
remove_on=0;
remove_threshold=40;

% thresholding
num_ca=1;
signal_per_ion=80;
noise=20;
threshold=signal_per_ion/log(1+signal_per_ion/noise);
threshold2=130;%sqrt(2)*signal_per_ion*1;
if num_ca==1
    threshold2=sqrt(2)*signal_per_ion*10;
end


%% IMPORT DATA
clear data
data=[];
% 1st: no. of repetition; 2nd: AOM freq; 3rd: LIF
for i=1:max(size(filenames))
    clear data1
    data1 = h5read(filenames(i).name,'/datasets/results');

    % get experiment parameters
    % readout time
    time = h5readatt(filenames(i).name,'/arguments/','time_sideband_readout_us');

    % sbc freq (only the first one)
    freq_sbc = h5readatt(filenames(i).name,'/arguments/','freq_sideband_cooling_mhz_pct_list');
    freq_sbc = str2double(extractBetween(freq_sbc, '{', ':'));
    ampl_quench_pct = h5readatt(filenames(i).name,'/arguments/','ampl_quench_pct');  %pct

    data=[data; data1'];
end

% remove data with consecutive 0s
if remove_on==1
    data=remove_zero(data,remove_threshold);
end


%% PROCESS DATA
% get unique frequencies
freq=unique(data(:,1));

% conduct binary thresholding and average to get probability
for i=1: max(size(freq))   
              order1=data(data(:,1)==freq(i) ,2)>threshold...
               &data(data(:,1)==freq(i),2)<threshold2;
           order2=data(data(:,1)==freq(i) ,2)>threshold2;
           LIF(i)=mean(order1+2*order2)/num_ca;   
           LIF_err(i)=std(order1+2*order2)/sqrt(max(size(order1)))/num_ca;
   %LIF(i)=mean(data(data(:,1)==freq(i),2));   
end

% convert frequency values from FTW to MHz
% note: FTW is Frequency Tuning Word, 0xFFFFFFFF = 1 GHz
freq2 = 2*freq*2.32830644e-7;

% separate RSB and BSB
freq_carrier_aom_mhz = mean(freq2);
freq_r = freq2(freq2 < freq_carrier_aom_mhz);
freq_b = freq2(freq2 > freq_carrier_aom_mhz);
LIF_r = LIF(freq2 < freq_carrier_aom_mhz);
LIF_r_err = LIF_err(freq2 < freq_carrier_aom_mhz);
LIF_b = LIF(freq2 > freq_carrier_aom_mhz);
LIF_b_err = LIF_err(freq2 > freq_carrier_aom_mhz);

%% PLOT DATA
figure(Fignum);

% plot RSB
subplot(1,2,1);hold on
errorbar(freq_r, LIF_r,LIF_r_err,'*-','DisplayName',[filenames(end).name(5:9) 'r'])
ylabel('S state population')
ylim([0.0, 1])
legend('Location','southwest')

% plot BSB
subplot(1,2,2);hold on;
errorbar(freq_b, LIF_b,LIF_b_err,'*-','DisplayName',[filenames(end).name(5:9) 'b'])
xlabel('Frequency -411041622@wavemeter (MHz)')
ylabel('S state population')
ylim([0.0, 1])
legend('Location','southwest')


%% FIT DATA
% guess fitting parameters
% guess RSB and BSB frequencies based on minima
minb_freq = mean(freq_b(find(LIF_b == min(LIF_b))) / 2);
minr_freq = mean(freq_r(find(LIF_r == min(LIF_r))) / 2);

% get carrier as mean of RSB and BSB freqs
guess_carrier = (minr_freq + minb_freq) / 2;

% guess offset as maximum-median
guess_bsb_offset = 1 - mean([max(LIF_b), median(LIF_b)]);
guess_rsb_offset = 1 - mean([max(LIF_r), median(LIF_r)]);

% guess ampl as minimum-offset
guess_bsb_ampl = (1 - min(LIF_b)) - guess_bsb_offset;
guess_bsb_ampl = 2 * asin(sqrt(guess_bsb_ampl)) / (pi * time);
guess_rsb_ampl = (1 - min(LIF_r)) - guess_rsb_offset;
guess_rsb_ampl = 2 * asin(sqrt(guess_rsb_ampl)) / (pi * time);

% fit BSB
fo_b = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, (2 * minb_freq) - 0.02, time, 0],...
               'Upper',[1, (2 * minb_freq) + 0.02, time, 0.3],...
               'StartPoint',[guess_bsb_ampl, 2 * minb_freq,  time, guess_bsb_offset]);
ftb = fittype('1 - a^2 / ((x-b)^2 + a^2) * (sin(2*pi*sqrt((x-b)^2 + a^2) * c/2))^2 - d','options',fo_b);
[curveb,gofb] = fit(freq_b,LIF_b',ftb);

% fit RSB
fo_r = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, (2 * minr_freq) - 0.02, time, 0],...
               'Upper',[1, (2 * minr_freq) + 0.02, time, 0.3],...
               'StartPoint',[guess_rsb_ampl, 2 * minr_freq,  time, guess_rsb_offset]);
ftr = fittype('1 - a^2 / ((x-b)^2 + a^2) * (sin(2*pi*sqrt((x-b)^2 + a^2) * c/2))^2 - d','options',fo_r);
[curver,gofr] = fit(freq_r,LIF_r',ftr);

% plot fitted RSB
figure(Fignum);hold on
subplot(1,2,1);plot(curver)
ylim([0.0,1]);ylabel('S state population')

% plot fitted BSB
subplot(1,2,2);plot(curveb)
xlabel('Frequency -411041622@wavemeter (MHz)')
ylabel('S state population');ylim([0.0,1]);legend('hide');


%% EXTRACT VALUES
% process RSB/BSB ratio to get phonon number
%R=(pi*curver.a)^2*curver.c^2/((pi*curveb.a)^2*curveb.c^2);
R = (1 - min(curver(freq_r)) - curver.d)./ (1 - min(curveb(freq_b)) - curveb.d);
n = R / (1 - R);

% state final phonon number on plot
title(sprintf('n = %.3f', round(100*n)/100));

% print relevant values
sprintf('RSB: %.4f \t BSB: %.4f\nQuench: %.2f%% \t SBC: %.4f%', ...
    curver.b / 2., curveb.b / 2., ampl_quench_pct, freq_sbc)
sprintf('n=%.3f', round(100. * n) / 100.);
