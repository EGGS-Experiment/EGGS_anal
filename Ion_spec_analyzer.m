%analysis quadrupole transition with 729 on
clear

filenames=[dir('*50329*.h5')  ];
Fignum=1095;

remove_on=0;
remove_threshold=100;

Ca_carrier_freq=206.124;%MHz
secular=1.088; %MHz

num_ca=1;
signal_per_ion=110;
noise=14;
threshold=55;signal_per_ion/log(1+signal_per_ion/noise);
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
    eggs_ht=h5readatt(filenames(i).name,'/arguments/','time_eggs_heating_ms');  %ms
    eggs_att=h5readatt(filenames(i).name,'/arguments/','att_eggs_heating_db');  %dB
end

if remove_on==1
    data=remove_zero(data,remove_threshold);
end

% todo: guess carrier freq
% todo: guess sideband freq


AOM_freq=unique(data(:,1));
carrier_freq_all=data(:,4);
secular_freq_all=data(:,5);

carrier_freq=unique(carrier_freq_all);
secular_freq=unique(secular_freq_all);

num_AOM=max(size(AOM_freq)) ;
num_carrier=max(size(unique(data(:,4))));
num_secular=max(size(unique(data(:,5))));


all_LIF_b=[];
all_LIF_r=[];

% 
% in this code, secular and carrier are flipped in order to quickly plot
% out data
    for m=1: num_secular  %(offset)
        clear LIF
        LIF=[];
        for i=1: num_AOM
            for k=1: num_carrier
                index=data(:,1)==AOM_freq(i) ...
                    &data(:,4)==carrier_freq(k)...
                    &data(:,5)==secular_freq(m);

               order1=data(index,2)>threshold &data(index ,2)<threshold2;
               order2=data(index ,2)>threshold2;
               LIF(i,k)=mean(order1+2*order2)/num_ca;   
               LIF_err(i,k)=std(order1+2*order2)/num_ca/sqrt(max(size(order1)));
            end
        end
        for k=1:num_carrier
           
            freq2=2*(AOM_freq)*2.32830644e-7;
            freq_r=freq2(freq2<Ca_carrier_freq);
            LIF_r=min(LIF(freq2<Ca_carrier_freq,k));
           % LIF_r_err= LIF_err(freq2<Ca_carrier_freq,k);

            freq_b=freq2(freq2>Ca_carrier_freq);
            LIF_b=min(LIF(freq2>Ca_carrier_freq,k));
           %LIF_b_err= LIF_err(freq2>Ca_carrier_freq,k);
           all_LIF_r(k)=1-LIF_r;
            figure(m+Fignum);subplot(1,2,1);hold on
            plot(freq_r,LIF(freq2<Ca_carrier_freq,k),'ro')

               xlabel('Sideband Detuning (MHz)')

              ylim([0 1])
            %legend 
            all_LIF_b(k)=1-LIF_b;
             figure(m+Fignum);subplot(1,2,2);hold on
             plot(freq_b,LIF(freq2>Ca_carrier_freq,k),'bo')
            ylim([0 1])

              xlabel('Sideband Detuning (MHz)')

             %legend 

        end
        offset(m)=secular_freq(m);
        ratio(m)=all_LIF_r(k)./all_LIF_b(k);


    end

figure;plot(offset,coherent_c(ratio),'o-')
