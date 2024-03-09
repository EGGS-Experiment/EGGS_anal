%analysis quadrupole transition with 729 on
clear

% filenames=[dir('*20069*.h5') dir('*20070*.h5') dir('*20071*.h5') dir('*20072*.h5')  dir('*20073*.h5')  ] ;
%filenames=[dir('*20076*.h5')   dir('*20077*.h5')   dir('*20078*.h5')  dir('*20079*.h5') ] ;
date_path = '\\eric.physics.ucla.edu\groups\motion\Data\2024-03\2024-03-05'; 
filenames=[dir(fullfile(date_path, '*25815*.h5'))]
Fignum=1094;

remove_on=0;
remove_threshold=100;

Ca_carrier_freq=206.816;%MHz
secular=1.101; %MHz

num_ca=1;
signal_per_ion=140;
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
    file = fullfile(date_path, filenames(i).name); 
    data1 = h5read(file,'/datasets/results');
    data=[data; data1'];
%     parameter_all=h5read(file,'/expid');
%     parameter=split(parameter_all,',');
    time=h5readatt(file,'/arguments/','time_eggs_heating_ms');
end

if remove_on==1
    data=remove_zero(data,remove_threshold);
end


AOM_freq=unique(data(:,1));
carrier_freq_all=data(:,4);
secular_freq_all=data(:,3);

carrier_freq=unique(carrier_freq_all);
secular_freq=unique(secular_freq_all);

num_AOM=max(size(AOM_freq)) ;
num_carrier=max(size(unique(data(:,4))));
num_secular=max(size(unique(data(:,3))));
% 


    for m=1: num_secular
        clear LIF
        LIF=[];
        for i=1: num_AOM
            for k=1: num_carrier
                index=data(:,1)==AOM_freq(i) ...
                    &data(:,4)==carrier_freq(k)...
                    &data(:,3)==secular_freq(m);

               order1=data(index,2)>threshold &data(index ,2)<threshold2;
               order2=data(index ,2)>threshold2;
               LIF(i,k)=mean(order1+2*order2)/num_ca;   
               LIF_err(i,k)=std(order1+2*order2)/num_ca/sqrt(max(size(order1)));
            end
        end
        for k=1:num_carrier
           
            freq2=2*(AOM_freq)*2.32830644e-7;
            freq_r=freq2(freq2<Ca_carrier_freq);
            LIF_r=LIF(freq2<Ca_carrier_freq,k);
            LIF_r_err= LIF_err(freq2<Ca_carrier_freq,k);


            freq_b=freq2(freq2>Ca_carrier_freq);
            LIF_b=LIF(freq2>Ca_carrier_freq,k);
           LIF_b_err= LIF_err(freq2>Ca_carrier_freq,k);
            
          figure(m+Fignum);subplot(2,1,1);hold on
            plot3(freq_r,1e-6*carrier_freq(k)*ones(1,max(size(freq_r))),1-LIF_r,'ro','DisplayName',num2str(carrier_freq(k)))  

           view(90,0)
             xlabel('Sideband Detuning (MHz)')
             ylabel('EGGs Frequency (MHz)')
            zlabel('D state')
            zlim([0 1])
            legend 

            figure(m+Fignum);subplot(2,1,2);hold on
            
            plot3(freq_b,1e-6*carrier_freq(k)*ones(1,max(size(freq_b))),1-LIF_b,'bo','DisplayName',num2str(carrier_freq(k)))

            view(90,0)
            xlabel('Sideband Detuning (MHz)')
            ylabel('EGGs Frequency (MHz)')

            zlabel('D state')
            zlim([0 1])
            legend 

        end
    end


% extract rabi freq
% extract bsb/rsb ratio
sb_r=1-LIF(1,:);
sb_b=1-LIF(2,:);
sb_ratio=(1-LIF(1,:))./(1-LIF(2,:));
tickle_freq_khz=(carrier_freq*1e-3);

% prepare fits: a, b, c, d
% a (rabi frequency): 
% b (center freq/resonance): x-value @ index of max y
% c (evolution time): time from exp submission
% d (offset): median val of points
param_b_start=tickle_freq_khz(sb_r==max(sb_r));
param_d_start=mean(sb_r);

fit_params = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, param_b_start-5, time, 0],...
               'Upper',[0.1, param_b_start+5, time, 1],...
               'StartPoint',[0.0005, param_b_start,  time, param_d_start]);
fit_params_full = fittype('1-a^2/((x-b)^2+a^2)*(sin(2*pi*sqrt((x-b)^2+a^2)*c/2))^2-d','options',fit_params);
[curveb,gofb] = fit(tickle_freq_khz, sb_r.', fit_params_full);

% plot resultant fits
%figure(1+Fignum);subplot(2,1,1);hold on
%plot(curveb)

