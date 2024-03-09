%% LINEWIDTH MEASUREMENT (397nm)

%% CONFIGURE
date_path = '\\eric.physics.ucla.edu\groups\motion\Data\2024-03\2024-03-08'; 
rid_dj = '51122';

rid_str = strcat('*', rid_dj, '*.h5');
filenames=[dir(fullfile(date_path, rid_str))];

fitfunc='g'; % g: gaussian fit; b: beta fit;
rf=19.057;
Fignum=1025;


%% PROCESS
% create data structures
clear data
data=[];

% import data
% 1st: no. of repetition; 2nd: AOM freq; 3rd: LIF
for i=1:max(size(filenames))
    file = fullfile(date_path, filenames(i).name); 
    clear data1
    data1 = h5read(file,'/datasets/results');
    data=[data; data1'];
end

% 
ion_on_all=data(data(:,2)==1,:);
ion_off_all=data(data(:,2)==0,:);       
LIF_on=ion_on_all(:,3);
LIF_off=ion_off_all(:,3);
RF=unique(data(:,1));


%% PROCESS DATA
freq=[];
LIF_mean=[];
LIF_err=[];
for i=1:max(size(unique(RF)))
    clear order
    %freq(i)=RF(i);
    freq(i)=RF(i)*2.38230644e-1;
    LIF_mean(i)=mean(LIF_on(ion_on_all(:,1)==RF(i))-LIF_off(ion_off_all(:,1)==RF(i)));
    %LIF_mean(i)=sum(LIF_on(ion_on_all(:,1)==RF(i)));
    LIF_err(i)=std(LIF_on(ion_on_all(:,1)==RF(i))-LIF_off(ion_off_all(:,1)==RF(i)))./sqrt(max(size(LIF_on)));
    %norm(i)=abs(mean(PD(order));
    
end

%% COLLATE DATA
% tmp remove
file1 = fullfile(date_path, filenames(1).name); 
data20 = h5read(file1,'/datasets/res_bgr');
data21 = h5read(file1,'/datasets/res_signal');
data3 = data21;

data3(1, :) = data3(1, :) * 2.;
data3(2, :) = data3(2, :) - data20(2, :);
data3 = data3.';

figure(Fignum);hold on;

last_file = fullfile(date_path, filenames(end).name); 
plot(data3(:,1), data3(:,2), 'o', 'DisplayName', [last_file(67:71)])
% tmp remove

figure(Fignum);hold on;
errorbar(2*freq*10^-6, LIF_mean, LIF_err, 'o', 'DisplayName',[last_file(67:71)])
xlabel('Abs. Freq. (MHz)');ylabel('counts/20us');legend('Location','northwest')

%% FITTING
xxxx=2*freq*10^-6;
yyyy=LIF_mean;
if fitfunc=='g'
    fo_b = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[0, 0, 0, 0],...
                   'Upper',[1, 1000, Inf, 1],...
                   'StartPoint',[0.1, 250, 30, 0.1]);
    ftb = fittype('a*exp(-1/2*(x-b)^2/c^2)+d','options',fo_b);
    [curveb,gofb] = fit(xxxx',yyyy',ftb);
    hold on;plot(curveb)
    sprintf('FWHM:\t\t%.2f MHz\nCenter:\t%.2f MHz', curveb.c*2.355, curveb.b)
elseif fitfunc=='b'
    P= @(center,beta,a,x)a*(besselj(-2,beta)^2./((x-2*rf-center).^2+10^2)+...
    besselj(-1,beta)^2./((x-rf-center).^2+10^2)+...
    besselj(0,beta)^2./((x-center).^2+10^2)+...
    besselj(1,beta)^2./((x+rf-center).^2+10^2)+...
     besselj(2,beta)^2./((x+rf*2-center).^2+10^2));
    fo_b = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[110, 0, 10],...
               'Upper',[270.1, 3, 300],...
               'StartPoint',[232, 1, 127]);
    ftb = fittype(P,'options',fo_b);
    [curveb,gofb] = fit(xxxx',yyyy',ftb);
    hold on;plot(curveb)
    sprintf('betal:\t\t%.2f MHz\nCenter:\t%.2f MHz', curveb.beta, curveb.center)
end
