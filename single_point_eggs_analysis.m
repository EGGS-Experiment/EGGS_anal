%analysis quadrupole transition with 729 on
%clear
%close (linspace(1097,1117,21))

filenames=[dir('*50329*.h5')  ];
Fignum=1098;

remove_on=0;
remove_threshold=100;

Ca_carrier_freq=205.6;%MHz
secular=0.770;%1.092; %MHz

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
    eggs_ht=h5readatt(filenames(i).name,'/arguments/','time_eggs_heating_ms');  %ms
    eggs_att=h5readatt(filenames(i).name,'/arguments/','att_eggs_heating_db');  %dB
    time=h5readatt(filenames(i).name,'/arguments/','time_sideband_readout_us');  %dB
    
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
            LIF_r=(LIF(freq2<Ca_carrier_freq,k));
            LIF_r_err= LIF_err(freq2<Ca_carrier_freq,k);

            freq_b=freq2(freq2>Ca_carrier_freq);
            LIF_b=(LIF(freq2>Ca_carrier_freq,k));
            LIF_b_err= LIF_err(freq2>Ca_carrier_freq,k);




            
           
            figure(k+Fignum);
            subplot(1,2,1);hold on
            plot(freq_r,LIF(freq2<Ca_carrier_freq,k),'or')
            ylim([0 1])
             xlabel('freq(MHz)')   
             ylabel('D state')
            
            figure(k+Fignum);subplot(1,2,2);hold on
            plot(freq_b,LIF(freq2>Ca_carrier_freq,k),'ob')
            xlabel('freq(MHz)')
            ylabel('D state')        
            ylim([0 1])    
            all_LIF_r(k)=1-LIF_r;
            all_LIF_b(k)=1-LIF_b;
    

        end
    end
%%%%%%%%%% FIT DATA %%%%%%%%%%
            % guess fitting parameters
            % guess RSB and BSB frequencies based on minima

        
figure(m+2*Fignum);
subplot(2,2,1);hold on
plot(carrier_freq*1e-6,all_LIF_r','or')

ylim([0 1])
 xlabel('secular freq(MHz)')   
 ylabel('D state')

figure(m+2*Fignum);subplot(2,2,2);hold on
plot(carrier_freq*1e-6,all_LIF_b','ob')
xlabel('secular freq(MHz)')
ylabel('D state')        
ylim([0 1])    
    
    
% fitting rsb with sinc^2

all_ratio=all_LIF_r./all_LIF_b;


all_ratio(all_ratio<0)=0;
all_ratio(all_ratio>0.791976610483122)=0.791976610483122;
all_phonon=coherent_c(all_ratio);
%all_phonon=squeezing_s(all_ratio);

figure(m+2*Fignum);subplot(2,2,3);plot(carrier_freq*1e-6,all_ratio,'o')
title('ratio RSB/BSB')
center=carrier_freq(all_ratio==max(all_ratio))*1e-6;
center=mean(center);
% fitting phonon with sinc^2
fo_phonon = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,center-0.01,1e3*eggs_ht 0],...
               'Upper',[2*pi*20e-3 center+0.01 1e3*eggs_ht 0.300],...
               'StartPoint',[2*pi*4e-3 center 1e3*eggs_ht 0.10]);
ftphonon = fittype('2*a^2/(x-b+1e-10)^2*(sin(2*pi*(x-b+1e-10)*c/2))^2+d','options',fo_phonon);
[curvephonon,gofphonon] = fit(carrier_freq*1e-6,all_phonon',ftphonon);

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



figure(m+2*Fignum)
title([filenames(end).name(5:9)])
subplot(2,2,4);hold on
plot(carrier_freq*1e-6,all_phonon,'o')
plot(carrier_freq*1e-6,curvephonon(carrier_freq*1e-6),'-')
title('phonon')
xlabel('secular freq(MHz)')
ylabel('Phonon') 

sprintf(['Rabi freq of heating is 2pi*',num2str(curveb.a*1e3),'kHz'])
%sprintf('Final phonon (1sigma): %s %s %s',num2str(curveb(center)-curveb.d),'\pm', num2str(err_phonon) )
sprintf('Final phonon (1sigma): %s %s %s',num2str(curvephonon(center)-curvephonon.d),'\pm', num2str(err_phonon) )

%(all_std(2,1)-all_std(1,1))/4  % 1sigma standard error

close (linspace(1097,1117,21));