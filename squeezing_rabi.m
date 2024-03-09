%% SQUEEZING RABI - FIT RABI FLOPPING FOR SQUEEZED/DISPLACED STATES
%figure;plot(data(:,2))
%figure;histogram(data(:,2))

%% CONFIGURE
date_path = '\\eric.physics.ucla.edu\groups\motion\Data\2024-03\2024-03-08'; 
rid_dj = '51120';

rid_str = strcat('*', rid_dj, '*.h5');
filenames=[dir(fullfile(date_path, rid_str))];


clear data data1 state delay LIF
Fignum=126;
sideband=0;
remove_on=0;
remove_threshold=50;

num_ca=1;
signal_per_ion=180;
noise=23.0;

delay_col = 1;   %  col:1 if run from experiment:rabi flopping ; col6 if run from qvsa; col7 if run from squeezing exp.
state_type = 't'; % 't' for thermal, 's' for squeezed, 'c' for coherent


%% PREPARE DATA
% set count thresholds
threshold=signal_per_ion/log(1+signal_per_ion/noise);
threshold2=170;sqrt(2)*signal_per_ion*1;
if num_ca==1
    threshold2=10*sqrt(2)*signal_per_ion*1;
end

% extract data
% 1st: scan time (s); 2nd: AOM freq; 3rd: Calibration; 4th: Fluorescence
data=[];
last_file = fullfile(date_path, filenames(end).name);
data1 = h5read(last_file,'/datasets/results');
% phase_rsb_turns = h5readatt(last_file,'/arguments/','phase_eggs_heating_rsb_turns_list');
% phase_bsb_turns = h5readatt(last_file,'/arguments/','phase_eggs_heating_bsb_turns');
data=[data; data1'];

% remove consistent zero counts
%figure;plot(data(:, 2))
if remove_on==1
    data=remove_zero(data,remove_threshold);
end

% plot only BSB points
% center_tmp = mean(data(:,1));
% data=data(data(:,1)>=center_tmp, :);


%% PROCESS DATA
% get rabi flopping times (i.e. x-axis) and photon counts (y-axis)
delay=unique(data(:,delay_col));
counts=data(:,2);

% threshold counts
for i=1: max(size(delay))
    LIF(i)=mean(data(data(:,delay_col)==delay(i) ,2));

    order1=data(data(:,delay_col) == delay(i) ,2)>threshold & data(data(:,delay_col) == delay(i) ,2)<threshold2; % threshold between zero and one ions
    order2=data(data(:,delay_col) == delay(i) ,2)>threshold2;% threshold between one and two ions
    state(i)=mean(order1+2*order2)/num_ca;
end    


%% PLOT RAW DATA
figure(Fignum);hold on;
last_file = fullfile(date_path, filenames(end).name);
plot(delay*1e-6, (1-state), 'DisplayName', [last_file(67:71)])
xlabel('Time (ms)')
ylabel('D state population')
legend


%% FIT DATA
% guess parameters
guess_rabi_wkhz = 115.0;
guess_detuning_wkhz = 0.0;
guess_ampl = 0.88;
guess_param = 0.2;
guess_time_decay_ms = 0.92;
start_param_arr = [2*pi*1e3 * guess_rabi_wkhz, guess_ampl, guess_detuning_wkhz * 2*pi*1e3, ...
    guess_param,  guess_time_decay_ms * 1e-3];
% start_param_arr = [2*pi*1e3 * 275, 0.95, 0, 3,  0.9e-3];

% fit!
fo_r = fitoptions('Method','NonlinearLeastSquares',...
                           'Lower',[2*pi*80e3, 0.6, 0., 0., 0.4e-3],...
                           'Upper',[2*pi*200e3, 1, 2*pi*6e3, 20., 2.5e-3],...
                           'StartPoint', start_param_arr);
ftr = fittype(@(Omega0, amp, detuning, nbar, t_dec, x) ...
    rabi_fitting_quantum_state(x, Omega0, amp, detuning, nbar, t_dec, state_type),'options',fo_r);
[curver,gofr] = fit(delay*1e-9, (1-state)', ftr);

% calculate errors
err_all=confint(curver);
err_n=abs(err_all(1,4)-err_all(2,4))/4;


%% FORMAT RESULTS
% get times and plot16
tt=(0: 0.001: max(delay*1e-6));
plot(tt, curver(tt * 1e-3))

if state_type=='t'
    state_str = 'nbar';
elseif state_type == 'c'
    state_str = '\alpha^2';
elseif state_type == 's'
    state_str = 'sq';
end

% format text
params_fit1 = sprintf('%s:  %.4f +/- %.4f', state_str, curver.nbar, err_n);
params_fit2 = sprintf('\\Omega_{0}:  %.2f kHz \n t_{dec}:  %.3f ms', curver.Omega0 / (2000.*pi), curver.t_dec * 1000.);
text(0.55, 0.15, params_fit1)
text(0.55, 0.05, params_fit2)
% sprintf('RSB: %.3f \t BSB: %.3f\nSum:\t%.3f', phase_rsb_turns, phase_bsb_turns, phase_rsb_turns+phase_bsb_turns)
