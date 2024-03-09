%% RABI FLOPPING - PROCESS
%figure;plot(data(:,2))
%figure;histogram(data(:,2))

%% CONFIGURE
date_path = '\\eric.physics.ucla.edu\groups\motion\Data\2024-03\2024-03-08'; 
rid_dj = '51120';

rid_str = strcat('*', rid_dj, '*.h5');
filenames=[dir(fullfile(date_path, rid_str))];

clear data state  delay LIF
data=[];

Fignum=124;
sideband=0;

remove_on=0;
remove_threshold=100;

num_ca=1;
signal_per_ion=160;
noise=25;

threshold=signal_per_ion/log(1+signal_per_ion/noise);
threshold2=200;sqrt(2)*signal_per_ion*1;

if num_ca==1
    threshold2=10*sqrt(2)*signal_per_ion*1;
end


%% IMPORT DATASET
% 1st: scan time (s); 2nd: AOM freq; 3rd: Calibration; 4th: Fluorescence
clear data1
last_file = fullfile(date_path, filenames(end).name); 
data1 = h5read(last_file,'/datasets/results');
ampl_probe_spinpol_pct = h5readatt(last_file,'/parameters','beams.ampl_pct.ampl_probe_spinpol_pct');  %pct
ampl_doppler_pct = h5readatt(last_file,'/parameters','beams.ampl_pct.ampl_pump_cooling_pct');  %pct
ampl_repump_cooling_pct = h5readatt(last_file,'/parameters','beams.ampl_pct.ampl_repump_cooling_pct');  %pct
ampl_repump_qubit_pct = h5readatt(last_file,'/parameters','beams.ampl_pct.ampl_repump_qubit_pct');  %pct
data=[data; data1'];


%% PROCESS DATASET
if remove_on==1
    data=remove_zero(data,remove_threshold);
end

delay=unique(data(:,1));
counts=data(:,2);


for i=1: max(size(delay))
    LIF(i)=mean(data(data(:,1)==delay(i) ,2));

    order1=data(data(:,1)==delay(i) ,2)>threshold&data(data(:,1)==delay(i) ,2)<threshold2; % threshold between zero and one ions
    order2=data(data(:,1)==delay(i) ,2)>threshold2;% threshold between one and two ions
    state(i)=mean(order1+2*order2)/num_ca;
end


%% PLOT FIGURE
%figure;plot(delay*1e3,1-state)
last_file = fullfile(date_path, filenames(end).name); 
figure(Fignum);hold on;plot(delay*1e-6,1-state,'DisplayName',[last_file(67:71)])
%figure(Fignum);hold on;plot(delay*1e3,1-state,'DisplayName',[filenames(end).name(5:9)])
xlabel('Time (ms)')
ylabel('D state population')
legend 

