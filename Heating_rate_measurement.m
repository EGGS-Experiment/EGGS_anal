% Heating_rate_measurement.m: Wait after SBC to measure the heating rate.
clear

% file/display
date_path = '\\eric.physics.ucla.edu\groups\motion\Data\2024-03\2024-03-05'; 
filenames=[dir(fullfile(date_path, '*50555*.h5'))]
Fignum=20932;

% remove data consisting of consecutive 0s
remove_on=0;
remove_threshold=100;

% thresholding
num_ca=1;
signal_per_ion=165;
noise=55;
threshold=signal_per_ion/log(1+signal_per_ion/noise);
threshold2=sqrt(2)*signal_per_ion*1;
if num_ca==1
    threshold2=sqrt(2)*signal_per_ion*10;
end


%%
%%%%%%%%%% IMPORT DATA %%%%%%%%%%
clear data
data=[];
% 1st: AOM freq; 2nd: LIF; 3rd: waiting time
for i=1:max(size(filenames))
    clear data1
    file = fullfile(date_path, filenames(i).name); 
    data1 = h5read(file,'/datasets/results');
    data=[data; data1'];
    data=data(data(:,1)~=0,:);

    % get experiment parameters
    % readout time
    time = h5readatt(file,'/arguments/','time_sideband_readout_us');
    time_readout = h5readatt(file,'/arguments/','time_sideband_readout_us');

    % sbc freq (only the first one)
    freq_sbc = h5readatt(file,'/arguments/','freq_sideband_cooling_mhz_pct_list');
    freq_sbc = str2double(extractBetween(freq_sbc, '{', ':'));
end

% remove data with consecutive 0s
if remove_on==1
    data=remove_zero(data,remove_threshold);
end


%%
%%%%%%%%%% PROCESS DATA %%%%%%%%%%
% get unique frequencies
freq = unique(data(:,1));
data(:,3) = data(:,3) * 1e-6; % change ns to ms
waiting_time_all = data(:,3);
waiting_time = unique(waiting_time_all);
num_waiting = max(size(unique(data(:,3))));

% guess carrier as average of maximum and minimum frequencies
guess_carrier_mhz = 0.5 * (min(freq) + max(freq)) * (2 * 2.32830644e-7);
guess_secular_mhz = max(freq) * (2 * 2.32830644e-7) - guess_carrier_mhz;


% prepare fitting for RSB/BSB
% fo_r = fitoptions('Method','NonlinearLeastSquares',...
%                'Lower',[0,      guess_carrier_mhz-secular,   time,   0],...
%                'Upper',[0.5,    guess_carrier_mhz-secular,  time,    0.2],...
%                'StartPoint',[0.0005, carrier_freq-secular, time, 0.1]);
fo_r = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, 0, time, 0],...
               'Upper',[1, Inf, time, 0.2],...
               'StartPoint',[0.0005, guess_carrier_mhz-guess_secular_mhz, time, 0.1]);
fo_b = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, 0, time, 0],...
               'Upper',[1, Inf, time, 0.2],...
               'StartPoint',[0.0005, guess_carrier_mhz+guess_secular_mhz, time, 0.1]);
ftr = fittype('1-a^2/((x-b)^2+a^2)*(sin(2*pi*sqrt((x-b)^2+a^2)*c/2))^2-d','options',fo_r);
ftb = fittype('1-a^2/((x-b)^2+a^2)*(sin(2*pi*sqrt((x-b)^2+a^2)*c/2))^2-d','options',fo_b);


%%
% PLOT & FIT HEATING RATES
for pp=1:num_waiting

    %% PREPARE VALUES FOR FITTING
    % extract values from dataset corresponding to given wait time
    for i=1: max(size(freq))   
        order1=data(data(:,1)==freq(i) &data(:,3)==waiting_time(pp),2)>threshold...
           &data(data(:,1)==freq(i)&data(:,3)==waiting_time(pp) ,2)<threshold2;
       order2=data(data(:,1)==freq(i)&data(:,3)==waiting_time(pp) ,2)>threshold2;
       LIF(i)=mean(order1+2*order2)/num_ca;   
    end

    % convert units to MHz
    freq2 = 2 * freq * 2.32830644e-7;

    % separate RSB from BSB
    freq_r = freq2(freq2 < guess_carrier_mhz);
    freq_b = freq2(freq2 > guess_carrier_mhz);
    LIF_r = LIF(freq2 < guess_carrier_mhz);
    LIF_b = LIF(freq2 > guess_carrier_mhz);

    % guess fitting parameters
    % guess RSB and BSB frequencies based on minima
    minb_freq = mean(freq_b(find(LIF_b == min(LIF_b))) / 2);
    minr_freq = mean(freq_r(find(LIF_r == min(LIF_r))) / 2);
    
    % guess offset as maximum-median
    guess_bsb_offset = 1 - mean([max(LIF_b), median(LIF_b)]);
    guess_rsb_offset = 1 - mean([max(LIF_r), median(LIF_r)]);
    
    % guess ampl as minimum-offset
    guess_bsb_ampl = (1 - min(LIF_b)) - guess_bsb_offset;
    guess_bsb_ampl = 2 * asin(sqrt(guess_bsb_ampl)) / (pi * time);
    guess_rsb_ampl = (1 - min(LIF_r)) - guess_rsb_offset;
    guess_rsb_ampl = 2 * asin(sqrt(guess_rsb_ampl)) / (pi * time);

    % fit RSB and BSB
    fo_r = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, (2 * minr_freq) - 5, time, 0],...
               'Upper',[0.01, (2 * minr_freq) + 5, time, 1],...
               'StartPoint',[guess_rsb_ampl, 2 * minr_freq,  time, guess_rsb_offset]);
    fo_b = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, (2 * minb_freq) - 5, time, 0],...
               'Upper',[0.01, (2 * minb_freq) + 5, time, 1],...
               'StartPoint',[guess_bsb_ampl, 2 * minb_freq,  time, guess_bsb_offset]);
    ftr = fittype('1-a^2/((x-b)^2+a^2)*(sin(2*pi*sqrt((x-b)^2+a^2)*c/2))^2-d','options',fo_r);
    ftb = fittype('1-a^2/((x-b)^2+a^2)*(sin(2*pi*sqrt((x-b)^2+a^2)*c/2))^2-d','options',fo_b);
    [curver,gofr] = fit(freq_r,LIF_r',ftr);
    [curveb,gofb] = fit(freq_b,LIF_b',ftb);


    %% PLOT!
    figure(Fignum); hold on
    % RSB subplot
    subplot(num_waiting+1,2,2*pp-1);
    last_file = fullfile(date_path, filenames(end).name); 
    plot(freq_r,LIF_r,'*-','DisplayName',[last_file(67:71) 'r'])
    hold on; plot(curver)
    ylabel('S state')
    ylim([0,1])
    legend('Location','southeast')
    title('Waiting Time (ms)', waiting_time(pp))

    % BSB subplot
    subplot(num_waiting+1,2,2*pp);
    plot(freq_b,LIF_b,'*-','DisplayName',[last_file(67:71) 'b'])
    hold on; plot(curveb)
    xlabel('Frequency -411041616@wavemeter (MHz)')
    ylabel('S state')
    ylim([0,1])
    legend('Location','southeast')
    
    %R=(pi*curver.a)^2*curver.c^2/((pi*curveb.a)^2*curveb.c^2);
    R=(1-min(curver(freq_r))-curver.d)./(1-min(curveb(freq_b))-curveb.d);
    
    % extract values for phonon number calculation
    wait(pp) = waiting_time(pp);
    phonon(pp) = R / (1 - R);

    % std for phonon number
    all_std = confint(curveb);
    bsb_a.err(pp) = (all_std(2,1)-all_std(1,1)) / 4;

    all_std = confint(curver);
    rsb_a.err(pp) = (all_std(2,1) - all_std(1,1)) / 4;

    photon_err(pp) = sqrt(rsb_a.err(pp)^2*(2*curver.a*(curveb.a)^2/((curveb.a)^2-(curver.a)^2)^2)^2)...
        +sqrt(bsb_a.err(pp)^2*(2*curveb.a*(curver.a)^2/((curveb.a)^2-(curver.a)^2)^2)^2);

end

% plot summary figure with error bars
figure(Fignum);subplot(num_waiting + 1, 2, [2 * num_waiting + 1, 2 * num_waiting + 2]);
errorbar(wait,phonon,photon_err,'o')
xlabel('time (ms)')
ylabel('phonon')

% print results to log
phonon
