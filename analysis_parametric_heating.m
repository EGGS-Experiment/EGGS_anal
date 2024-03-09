%analysis quadrupole transition with 729 on
%clear

filenames=[dir('*36311*.h5')  ];
Fignum=1096;

remove_on=0;
remove_threshold=100;

Ca_carrier_freq=206.045;%MHz
secular=0.770*2;%1.092; %MHz

num_ca=1;
signal_per_ion=110;
noise=14;
threshold=65;signal_per_ion/log(1+signal_per_ion/noise);
threshold2=180;sqrt(2)*signal_per_ion;
if num_ca==1
    threshold2=sqrt(2)*signal_per_ion*10;
end

clear data
data=[];
% 1st: AOM freq; 2nd: LIF; ; 3th: carrier freq; 4th: secular frequency

for i=1:max(size(filenames))
    clear data1
    data1 = h5read(filenames(i).name,'/datasets/results');
    data=[data; data1'];
    %eggs_ht=h5readatt(filenames(i).name,'/arguments/','time_eggs_heating_ms');  %ms
    %eggs_att=h5readatt(filenames(i).name,'/arguments/','att_eggs_heating_db');  %dB
   % time=h5readatt(filenames(i).name,'/arguments/','time_sideband_readout_us');  %dB
    
end

if remove_on==1
    data=remove_zero(data,remove_threshold);
end

% todo: guess carrier freq
% todo: guess sideband freq


AOM_freq=unique(data(:,1));
carrier_freq_all=data(:,4);
secular_freq_all=data(:,3);

carrier_freq=unique(carrier_freq_all);
secular_freq=unique(secular_freq_all);

num_AOM=max(size(AOM_freq)) ;
num_carrier=max(size(unique(data(:,4))));
num_secular=max(size(unique(data(:,3))));
all_LIF_b=[];
all_LIF_r=[];

% 
% in this code, secular and carrier are flipped in order to quickly plot
% out data
    for m=1: num_secular
        clear LIF
        LIF=[];
        for i=1: num_AOM
            for p=1: num_secular
                index=data(:,1)==AOM_freq(i) ...
                    &data(:,3)==secular_freq(p);

               order1=data(index,2)>threshold &data(index ,2)<threshold2;
               order2=data(index ,2)>threshold2;
               LIF(i,p)=mean(order1+2*order2)/num_ca;   
               LIF_err(i,p)=std(order1+2*order2)/num_ca/sqrt(max(size(order1)));
            end
        end
        
           
            freq2=2*(AOM_freq)*2.32830644e-7;
            freq_r=freq2(freq2<Ca_carrier_freq);
            LIF_r=(LIF(freq2<Ca_carrier_freq,m));

            freq_b=freq2(freq2>Ca_carrier_freq);
            LIF_b=(LIF(freq2>Ca_carrier_freq,m));
 
        all_LIF_r(m)=1-LIF_r;
        all_LIF_b(m)=1-LIF_b;
    end
%%%%%%%%%% FIT DATA %%%%%%%%%%
            % guess fitting parameters
            % guess RSB and BSB frequencies based on minima

    
    
% fitting rsb with sinc^2

all_ratio=all_LIF_r./all_LIF_b;



figure(m+2*Fignum);subplot(2,2,3);plot(secular_freq*2.32830644e-7,all_ratio,'o')
title('ratio RSB/BSB')

