%% ISA FITTING
%analysis quadrupole transition with 729 on
clear

filenames=[dir('*47151*.h5')];
Fignum=1097;

remove_on=1;
remove_threshold=50;

Ca_carrier_freq=205.1570;%MHz
secular=0.759; %MHz

num_ca=1;
signal_per_ion=130;
noise=18;
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
    eggs_ht=2;%h5readatt(filenames(i).name,'/arguments/','time_eggs_heating_ms');  %ms
    eggs_att=26;%h5readatt(filenames(i).name,'/arguments/','att_eggs_heating_db');  %dB
    time=95;%h5readatt(filenames(i).name,'/arguments/','time_sideband_readout_us');  %us
    
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
            LIF_r=(LIF(freq2<Ca_carrier_freq,k));
           % LIF_r_err= LIF_err(freq2<Ca_carrier_freq,k);

            freq_b=freq2(freq2>Ca_carrier_freq);
            LIF_b=(LIF(freq2>Ca_carrier_freq,k));
           %LIF_b_err= LIF_err(freq2>Ca_carrier_freq,k);
           all_LIF_r(k)=mean(1-LIF_r);


                        % guess RSB and BSB frequencies based on minima
            minb_freq = mean(freq_b(find(LIF_b == min(LIF_b))) / 2);
            minr_freq = mean(freq_r(find(LIF_r == min(LIF_r))) / 2);
            


            % fit BSB
            fo_b = fitoptions('Method','NonlinearLeastSquares',...
                           'Lower',[0, (2 * minb_freq) - 4e-3, time-15, 0],...
                           'Upper',[10e-3, (2 * minb_freq) + 4e-3, time+15, 0.2],...
                           'StartPoint',[5e-3, 2 * minb_freq,  time, 0.05]);
            ftb = fittype('1 - a^2 / ((x-b)^2 + a^2) * (sin(2*pi*sqrt((x-b)^2 + a^2) * c/2))^2 - d','options',fo_b);
            [curveb,gofb] = fit(freq_b,LIF(freq2>Ca_carrier_freq,k),ftb);
            
            % fit RSB
            fo_r = fitoptions('Method','NonlinearLeastSquares',...
                           'Lower',[0, (2 * minr_freq) - 4e-3,curveb.c, 0],...
                           'Upper',[5e-3, (2 * minr_freq) + 4e-3, curveb.c, 0.2],...
                           'StartPoint',[1e-3, 2 * minr_freq,  curveb.c, curveb.d]);
            ftr = fittype('1 - a^2 / ((x-b)^2 + a^2)*(sin(2*pi*sqrt((x-b)^2 + a^2) * c/2))^2 - d','options',fo_r);
            [curver,gofr] = fit(freq_r,LIF(freq2<Ca_carrier_freq,k),ftr);


            figure(m+Fignum);subplot(1,2,1);hold on
            plot(freq_r,LIF(freq2<Ca_carrier_freq,k),'ro')
            plot(curver)


            xlabel('Sideband Detuning (MHz)')

            ylim([0 1])
            %legend 
            all_LIF_b(k)=mean(1-LIF_b);
            figure(m+Fignum);subplot(1,2,2);hold on
            plot(freq_b,LIF(freq2>Ca_carrier_freq,k),'bo')
            plot(curveb)
            ylim([0 1])

            xlabel('Sideband Detuning (MHz)')

             %legend 

        end
        offset(m)=secular_freq(m);
        ratio(m)=(1-curver(curver.b))/(1-curveb(curveb.b));


    end

figure;plot(offset,coherent_c(ratio),'o-')
