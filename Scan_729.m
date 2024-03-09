%analysis quadrupole transition with 729 on
%figure;plot(data(:,2))
%figure;histogram(data(:,2))
clear

%% CONFIGURE
date_path = '\\eric.physics.ucla.edu\groups\motion\Data\2024-03\2024-03-08'; 
rid_dj = '51127';

rid_str = strcat('*', rid_dj, '*.h5');
filenames=[dir(fullfile(date_path, rid_str))];

num_ca=1;
signal_per_ion=140;
noise=25;
Fignum=812;
plot_rf=0;

remove_on=0;
remove_threshold=100;


%% THRESHOLDING
threshold=signal_per_ion/log(1+signal_per_ion/noise);
%threshold2=170;
threshold2=sqrt(2)*signal_per_ion*1;
if num_ca==1
    threshold2=sqrt(2)*signal_per_ion*10;
end
        
%threshold=55;
%threshold2=140;


clear data
data=[];
data_rf=[];
% 1st: no. of repetition; 2nd: AOM freq; 3rd: LIF
for i=1:max(size(filenames))
    clear data1
    clear data_rf1
    file = fullfile(date_path, filenames(i).name); 
    data1 = h5read(file,'/datasets/results');
    %data_rf1 = h5read(file,'/datasets/adc_values');
    data=[data; data1'];
    %data_rf=[data_rf; data_rf1'];
end

% figure;plot(data(:,2))
if remove_on==1
    data=remove_zero(data,remove_threshold);
end
freq=unique(data(:,1));

for i=1: max(size(freq))
    
              order1=data(data(:,1)==freq(i) ,2)>threshold...
               &data(data(:,1)==freq(i),2)<threshold2;
           order2=data(data(:,1)==freq(i) ,2)>threshold2;
           LIF(i)=mean(order1+2*order2)/num_ca;   
   %LIF(i)=mean(data(data(:,1)==freq(i),2));
    
end

freq2=2*freq*2.32830644e-7;
figure(Fignum);hold on;plot(freq2,LIF,'DisplayName',[filenames(end).name(4:9)])
xlabel('Frequency -411041622@wavemeter (MHz)')
%legend(filenames(1).name(5:9), 'Location','southwest')
legend('Location','southwest')

% tmp remove
freq2(find(LIF==min(LIF))) / 2

