%analysis quadrupole transition with 729 on
clear

% filenames=[dir('*20069*.h5') dir('*20070*.h5') dir('*20071*.h5') dir('*20072*.h5')  dir('*20073*.h5')  ] ;
%filenames=[dir('*20076*.h5')   dir('*20077*.h5')   dir('*20078*.h5')  dir('*20079*.h5') ] ;
filenames=[dir('*24986*.h5')  ] ;
Fignum=100;

remove_on=1;
remove_threshold=100;

Ca_carrier_freq=206.822;%MHz
secular=1.101; %MHz

num_ca=1;
signal_per_ion=140;
noise=14;
threshold=signal_per_ion/log(1+signal_per_ion/noise);
threshold2=sqrt(2)*signal_per_ion;
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
%     parameter_all=h5read(filenames(i).name,'/expid');
%     parameter=split(parameter_all,',');
end

if remove_on==1
    data=remove_zero(data,remove_threshold);
end


AOM_freq=unique(data(:,1));
carrier_freq_all=data(:,3);
secular_freq_all=data(:,4);

carrier_freq=unique(carrier_freq_all);
secular_freq=unique(secular_freq_all);

num_AOM=max(size(AOM_freq)) ;
num_carrier=max(size(unique(data(:,3))));
num_secular=max(size(unique(data(:,4))));
% 
meow=[];

    for m=1: num_secular
        clear LIF
        LIF=[];
        for i=1: num_AOM
            for k=1: num_carrier
                index=data(:,1)==AOM_freq(i) ...
                    &data(:,3)==carrier_freq(k)...
                    &data(:,4)==secular_freq(m);

               order1=data(index,2)>threshold &data(index ,2)<threshold2;
               order2=data(index ,2)>threshold2;
               LIF(i,k)=mean(order1+2*order2)/num_ca;   
               LIF_err(i,k)=std(order1+2*order2)/num_ca/sqrt(max(size(order1)));
            end
        end
        for k=1:num_carrier
            %freq2=2*(AOM_freq);
            freq2=2*AOM_freq*2.32830644e-7;
            freq_r=freq2(freq2<Ca_carrier_freq);
            LIF_r=LIF(freq2<Ca_carrier_freq,k);
            LIF_r_err= LIF_err(freq2<Ca_carrier_freq,k);


            freq_b=freq2(freq2>Ca_carrier_freq);
            LIF_b=LIF(freq2>Ca_carrier_freq,k);
           LIF_b_err= LIF_err(freq2>Ca_carrier_freq,k);
            
            figure(m+Fignum);subplot(3,1,1);hold on
           % plot3(freq_r,carrier_freq(k)*ones(1,max(size(freq_r))),LIF_r,'ro')  
            errorbar(carrier_freq(k),1-LIF_r,LIF_r_err,'ro')  
            LIF_r_all(k)=1-LIF_r;
         %    view(0,0)
             xlabel('Sideband Detuning (MHz)')
             ylabel('EGGs Frequency (MHz)')
            zlabel('D state')
            zlim([0 1])
            title('secular frequency(MHz)', secular_freq(m))

            figure(m+Fignum);subplot(3,1,2);hold on
            
            %plot3(freq_b,carrier_freq(k)*ones(1,max(size(freq_b))),LIF_b,'bo')
            errorbar(carrier_freq(k),1-LIF_b,LIF_b_err,'bo')
            LIF_b_all(k)=1-LIF_b;
           % view(0,0)
            xlabel('Sideband Detuning (MHz)')
            ylabel('EGGs Frequency (MHz)')

            zlabel('D state')
            zlim([0 1])
            title('secular frequency(MHz)', secular_freq(m))
        end
        mean_b=0;
        mean_r=0;
        r_amp=(LIF_r_all-mean_r);
        b_amp=(LIF_b_all-mean_r);

        figure(m+Fignum);subplot(3,1,3);plot(carrier_freq,r_amp./b_amp...
        ./(1-r_amp./b_amp),'-o')


    end




% mean_b=0;%mean(LIF_b_all);
% mean_r=0;%mean(LIF_r_all);
%figure;plot(carrier_freq,LIF_r_all./LIF_b_all./(1-LIF_r_all./LIF_b_all))
% r_amp=(LIF_r_all-mean_r);
% b_amp=(LIF_b_all-mean_r);
%phonon_err= sqrt(LIF_r_err.^2.*(2.*r_amp.*(b_amp).^2./((b_amp).^2-(r_amp).^2).^2).^2)...
%           +sqrt(LIF_b_err.^2.*(2.*b_amp.*(r_amp).^2./((b_amp).^2-(r_amp).^2).^2).^2);

% 
% figure;errorbar(carrier_freq,r_amp./b_amp...
%     ./(1-r_amp./b_amp),phonon_err,'o')

figure;plot(carrier_freq,r_amp./b_amp...
    ./(1-r_amp./b_amp),'-o')
