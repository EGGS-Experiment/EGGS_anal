%% SQUEEZING RABI - FIT RABI FLOPPING FOR SQUEEZED/DISPLACED STATES
date_path = '\\eric.physics.ucla.edu\groups\motion\Data\2024-03\2024-03-07'; 
filenames=[dir(fullfile(date_path, '*47261*.h5'))];

clear data data1 state delay LIF
Fignum=122;
sideband=0;
remove_on=0;
remove_threshold=100;

num_ca=1;
signal_per_ion=100;
noise=15;

delay_col = 5;
state_type = 'c';

%% PROCESS DATA
% set count thresholds
threshold=58;signal_per_ion/log(1+signal_per_ion/noise);
threshold2=170;sqrt(2)*signal_per_ion*1;
if num_ca==1
    threshold2=10*sqrt(2)*signal_per_ion*1;
end

% extract data
% 1st: scan time (s); 2nd: AOM freq; 3rd: Calibration; 4th: Fluorescence
data=[];
last_file = fullfile(date_path, filenames(end).name);
data1 = h5read(last_file,'/datasets/results');
data=[data; data1'];
if remove_on==1
    data=remove_zero(data,remove_threshold);
end

% get rabi flopping times (i.e. x-axis)
delay=unique(data(:,delay_col));
% get photon counts
counts=data(:,2);

% threshold counts
for i=1: max(size(delay))
    LIF(i)=mean(data(data(:,delay_col)==delay(i) ,2));

    order1=data(data(:,delay_col) == delay(i) ,2)>threshold&data(data(:,delay_col) == delay(i) ,2)<threshold2; % threshold between zero and one ions
    order2=data(data(:,delay_col) == delay(i) ,2)>threshold2;% threshold between one and two ions
    state(i)=mean(order1+2*order2)/num_ca;
    state_err(i)=std(order1+2*order2)/sqrt(max(size(order1)))/num_ca;
end    


%% PLOT RAW DATA
%figure;plot(delay*1e3,LIF)
%figure;plot(delay*1e3,1-state)
figure(Fignum);hold on;
last_file = fullfile(date_path, filenames(end).name);
errorbar(delay*1e-6, (1-state),state_err,'o', 'DisplayName', [last_file(67:71)])
%figure(Fignum);hold on;plot(delay*1e3,1-state,'DisplayName',[last_file(67:71)])
xlabel('time (ms)')
xlabel('Time (ms)')
ylabel('D state population')
legend


%% FIT DATA
fo_r = fitoptions('Method','NonlinearLeastSquares',...
                           'Lower',[0, 0.98, 0 ,0.78, 0.4e-3],...
                           'Upper',[2*pi*475e3, 1.2, 2*pi*1e3, 1.75, 1.9e-3],...
                           'StartPoint',[2*pi*1e3 * 255, 0.95, 0, .93,  0.9e-3]);

ftr = fittype(@(Omega0,amp,detuning,nbar,t_dec,x) ...
    rabi_fitting_quantum_state(x,Omega0,amp,detuning,nbar, t_dec, state_type),'options',fo_r);
[curver,gofr] = fit(delay*1e-9,(1-state)',ftr);
tt=(0:0.001:1);
plot(tt,curver(tt*1e-3),'LineWidth',2)


%% FORMAT RESULTS
if state_type=='t'
    state_str = 'nbar';
elseif state_type == 'c'
    state_str = '\alpha';
elseif state_type == 's'
    state_str = 'sq';
end
err_all=confint(curver);
err_n=abs(err_all(1,4)-err_all(2,4))/4;
params_fit1 = sprintf('%s:  %.4f +/- %.4f', state_str, curver.nbar, err_n);
params_fit2 = sprintf('Omega0:  %.2f kHz \n t_{dec}:  %.3f ms', curver.Omega0 / (2000.*pi), curver.t_dec * 1000.);
text(0.55, 0.15, params_fit1)
text(0.55, 0.05, params_fit2)

first_file = fullfile(date_path, filenames(1).name);

sprintf(num2str(h5readatt(first_file,'/arguments/','phase_eggs_heating_rsb_turns')))
