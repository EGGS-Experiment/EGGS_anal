%Use this code when you run either EGGs Heating or Ion Spectrum Analyzer
%with multiple sideband values, a single readout time, multiple secular frequency values and one
%detuning i.e. on resonance.
%figure;plot(data(:,2))
%figure;histogram(data(:,2))
data=[];

%% CONFIGURE
date_path = '\\eric.physics.ucla.edu\groups\motion\Data\2024-03\2024-03-08'; 
rid_dj = '51120';

rid_str = strcat('*', rid_dj, '*.h5');
filenames=[dir(fullfile(date_path, rid_str))];
Fignum=1096;

% REMOVE ZERO
remove_on=0;
remove_threshold=50;

% FLUORESCENCE
num_ca=1;
noise=25;
signal_per_ion=140;

threshold = signal_per_ion/log(1+signal_per_ion/noise);
threshold2 = sqrt(2)*signal_per_ion;
if num_ca==1
    threshold2=sqrt(2)*signal_per_ion*10;
end


%% IMPORT DATASET
for i=1:max(size(filenames))
    clear data1
    file = fullfile(date_path, filenames(i).name); 
    data1 = h5read(file,'/datasets/results');
    data=[data; data1'];
    eggs_ht=h5readatt(file,'/arguments/','time_eggs_heating_ms');  %ms
    eggs_att=h5readatt(file,'/arguments/','att_eggs_heating_db');  %dB
    time=h5readatt(file,'/arguments/','time_sideband_readout_us');  %dB

    eggs_freq_carrier_mhz = h5readatt(file,'/arguments/','freq_eggs_heating_carrier_mhz_list');
    eggs_ampl_carrier_pct = h5readatt(file,'/arguments/','ampl_eggs_dynamical_decoupling_pct');
    eggs_ampl_rsb_pct = h5readatt(file,'/arguments/','ampl_eggs_heating_rsb_pct');
    eggs_ampl_bsb_pct = h5readatt(file,'/arguments/','ampl_eggs_heating_bsb_pct');
end

if remove_on==1
    data=remove_zero(data,remove_threshold);
end


%% PRE-PROCESS DATASET
aom_freq_index =        1;
carrier_freq_index =    4;
secular_freq_index =    3;
AOM_freq=unique(data(:,aom_freq_index));
carrier_freq_all=data(:,carrier_freq_index);
secular_freq_all=data(:,secular_freq_index);

Ca_carrier_freq = mean(AOM_freq)*2.e3/(2.^32.-1.);
carrier_freq=unique(carrier_freq_all);
secular_freq=unique(secular_freq_all);

num_AOM=max(size(AOM_freq)) ;
num_carrier=max(size(unique(data(:,carrier_freq_index))));
num_secular=max(size(unique(data(:,secular_freq_index))));
all_LIF_b=[];
all_LIF_r=[];

fit_freq_r_list=[];
fit_freq_b_list=[];


%% PLOT & FIT INDIVIDUAL POINTS
% 1st: AOM freq; 2nd: LIF; ; 3th: carrier freq; 4th: secular frequency
% note: secular and carrier values are flipped 
    for m=1: num_secular
        clear LIF
        LIF=[];

        % threshold & sort readout points
        for i=1: num_AOM
            for k=1: num_carrier
                index=data(:,aom_freq_index)==AOM_freq(i) ...
                    &data(:,carrier_freq_index)==carrier_freq(k)...
                    &data(:,secular_freq_index)==secular_freq(m);

               order1=data(index,2)>threshold & data(index ,2)<threshold2;
               order2=data(index ,2)>threshold2;
               LIF(i,k)=mean(order1+2*order2)/num_ca;   
               LIF_err(i,k)=std(order1+2*order2)/num_ca/sqrt(max(size(order1)));
            end
        end

        % fit readout spectrum at each sweep point
        for k=1:num_carrier
           
            % process readout spectrum points
            freq2=2*AOM_freq*2.32830644e-7;
            freq_r=freq2(freq2<Ca_carrier_freq);
            LIF_r=LIF(freq2<Ca_carrier_freq,k);
            LIF_r_err= LIF_err(freq2<Ca_carrier_freq,k);

            freq_b=freq2(freq2>Ca_carrier_freq);
            LIF_b=(LIF(freq2>Ca_carrier_freq,k));
            LIF_b_err= LIF_err(freq2>Ca_carrier_freq,k);

            minb_freq = mean(freq_b(find(LIF_b == min(LIF_b))) / 2);
            minr_freq = mean(freq_r(find(LIF_r == min(LIF_r))) / 2);

            % guess carrier as mean of RSB and BSB freqs
            guess_carrier = (minr_freq + minb_freq) / 2;
            guess_secular_mhz = 2*(minb_freq-guess_carrier);

            % guess offset as maximum-median
%             guess_bsb_offset = 1 - mean([min(LIF_b), median(LIF_b)]);
%             guess_rsb_offset = 1 - mean([min(LIF_r), median(LIF_r)]);
            guess_bsb_offset = 1 - max(LIF_b);
            guess_rsb_offset = 1 - max(LIF_r);

            % guess ampl as minimum-offset
            guess_bsb_ampl = (1 - min(LIF_b)) - guess_bsb_offset;
            guess_bsb_ampl = 2 * asin(sqrt(guess_bsb_ampl)) / (pi * time);
            guess_rsb_ampl = (1 - min(LIF_r)) - guess_rsb_offset;
            guess_rsb_ampl = 2 * asin(sqrt(guess_rsb_ampl)) / (pi * time);

            % fit BSB
            fo_b = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, (2*minb_freq) - 20e-3, time, 0],...
               'Upper',[2*pi*10e-3, (2*minb_freq) + 20e-3, time, 0.05],...
               'StartPoint',[guess_bsb_ampl, 2 * minb_freq,  time, guess_bsb_offset]);
            ftb = fittype('1 - a^2 / ((x-b)^2 + a^2) * (sin(2*pi*sqrt((x-b)^2 + a^2) * c/2))^2 - d','options',fo_b);
            [curveb,gofb] = fit(freq_b,LIF(freq2>Ca_carrier_freq,k),ftb);
            
            % fit RSB
            fo_r = fitoptions('Method','NonlinearLeastSquares',...
                           'Lower',[0, (2 * minr_freq) - 2e-3, time, 0],...
                           'Upper',[0.01, (2 * minr_freq) + 2e-3, time, 0.1],...
                           'StartPoint',[guess_rsb_ampl, 2 * minr_freq,  time, guess_rsb_offset]);
            ftr = fittype('1 - a^2 / ((x-b)^2 + a^2)*(sin(2*pi*sqrt((x-b)^2 + a^2) * c/2))^2 - d','options',fo_r);
            [curver,gofr] = fit(freq_r,LIF(freq2<Ca_carrier_freq,k),ftr);

            % plot RSB
            figure(k+Fignum);subplot(1,2,1);hold on
            plot(freq_r,LIF(freq2<Ca_carrier_freq,k),'or')
            plot(curver,'r');legend('hide')
            xlabel('freq(MHz)');ylabel('D state');ylim([0 1]);
            % plot BSB
            figure(k+Fignum);subplot(1,2,2);hold on
            plot(freq_b,LIF(freq2>Ca_carrier_freq,k),'ob')
            plot(curveb,'b');legend('hide');
            xlabel('freq(MHz)');ylabel('D state');ylim([0 1]);

            % extract amplitude from readout spectrum
            all_LIF_r(k)=1-curver(curver.b)-curver.d;
            all_LIF_b(k)=1-curveb(curveb.b)-curveb.d;
%             all_LIF_r(k)=curver.b;
%             all_LIF_b(k)=curveb.b;
            fit_freq_r_list(k) = curver.b;
            fit_freq_b_list(k) = curveb.b;
        end
    end

    
%% SUMMARY FIGURE FIT
% % guess fitting parameters
% % guess RSB and BSB frequencies based on minima
% 
% % guess secular
% guess_secular_mhz = 2*(minb_freq-guess_carrier);
% 
% % guess offset as maximum-median
% guess_bsb_offset = 1 - mean([max(LIF_b), median(LIF_b)]);
% guess_rsb_offset = 1 - mean([max(LIF_r), median(LIF_r)]);
% 
% % guess ampl as minimum-offset
% guess_bsb_ampl = (1 - min(LIF_b)) - guess_bsb_offset;
% guess_bsb_ampl = 2 * asin(sqrt(guess_bsb_ampl)) / (pi * time);
% guess_rsb_ampl = (1 - min(LIF_r)) - guess_rsb_offset;
% guess_rsb_ampl = 2 * asin(sqrt(guess_rsb_ampl)) / (pi * time);

% fit RSB summary
% guess_secular_eggs_mhz = ;
% secular = 0.7625;
secular = abs(0.5*(mean(fit_freq_r_list)-mean(fit_freq_b_list)));
fo_red = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, secular-0.025, 0.0001, 0.],...
               'Upper',[1, secular+0.025,  0.3, 0.3],...
               'StartPoint',[0.13, secular,  0.002, 0.12]);
ft_red = fittype('a*exp(-(x-b)^2/c^2)+d','options',fo_red);
[curve_red,gof_red] = fit(carrier_freq*1e-6, all_LIF_r',ft_red);
bg=curve_red.d;
center=curve_red.b;

% CONVERT RSB/BSB TO DISPLACEMENT
all_ratio=all_LIF_r./all_LIF_b;
all_ratio(all_ratio<0)=0;
all_ratio(all_ratio>0.791976610483122)=0.791976610483122;
all_phonon=coherent_c(all_ratio);
%all_phonon=squeezing_s(all_ratio);

% FIT <n> WITH SINC^2
fo_phonon = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, curve_red.b-0.01, 1e3*eggs_ht, 0],...
               'Upper',[2*pi*20e-3, curve_red.b+0.01, 1e3*eggs_ht, 0.300],...
               'StartPoint',[2*pi*4e-3, curve_red.b, 1e3*eggs_ht, 0.10]);
ftphonon = fittype('2*a^2/(x-b+1e-10)^2*(sin(2*pi*(x-b+1e-10)*c/2))^2+d','options',fo_phonon);
[curvephonon,gofphonon] = fit(carrier_freq*1e-6,all_phonon',ftphonon);

% PLOT RSB
figure(m+2*Fignum);subplot(2,2,1);hold on
plot(carrier_freq*1e-6,all_LIF_r','or')
plot(curve_red,'r');legend('hide')
xlabel('Secular Freq. (MHz)');ylabel('D state');%ylim([0, 1]);
% PLOT BSB
figure(m+2*Fignum);subplot(2,2,2);hold on
plot(carrier_freq*1e-6,all_LIF_b','ob')
xlabel('Secular Freq. (MHz)');ylabel('D state');%ylim([0, 1]);
% PLOT RSB/BSB
figure(m+2*Fignum);subplot(2,2,3);
plot(carrier_freq*1e-6,all_ratio,'o')
title('RSB/BSB')
% PLOT FITTED <n>
figure(m+2*Fignum);subplot(2,2,4);hold on
plot(carrier_freq*1e-6,all_phonon,'o')
plot(curvephonon)
legend('hide');title('<n>')
xlabel('secular freq(MHz)');ylabel('Phonon');


%% CALCULATE ERROR
all_std=confint(curvephonon);
curvephonon.a*1e3;  %  *2pi*kHz
curvephonon(center)-curvephonon.d;
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
    err_d=(all_std(2,4)-all_std(1,4))/4 ;
end
err_phonon=sqrt((4*pi^2*curvephonon.a*curvephonon.c^2*err_a)^2+(4*pi^2*curvephonon.a^2*curvephonon.c*err_c)^2+err_d^2);


%% PRINTOUT
% sprintf('Heating Rabi Freq:\t %.3f x2pi kHz\nFinal phonon (1sigma):\t %.3f +/- %.3f', ...
%     curveb.a*1e3, curvephonon(curvephonon.b)-curvephonon.d, err_phonon)

figure(m+3*Fignum);
histogram(fit_freq_r_list/2);
figure((m+1)+3*Fignum);
histogram(fit_freq_b_list/2);
mean(fit_freq_b_list)/2
mean(fit_freq_r_list)/2