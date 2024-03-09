freq_dipole=[50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 135., 140., 145., 150., 155., 165.];
%relative to 20dB
ampl_c=[6., 5.5, 4.25, 9.5, 30., 22.5, 14., 48., 37.5, 30., 29., 36., 77., 180., 440., 869.446, 370., 360., 800., 2200., 5130.5, 7915.62, 9391.25];
ampl_sb=40;

phonon=[0.9272, 1.0898, 1.0011, 0.8882, 0.809, 0.897, 0.8922, 0.7851, 1.0806, 0.7088, 0.7916, 0.7734, 0.7921, 0.6937, 0.9101, 0.8815, 1.1053, 0.9989, 0.7446, 0.8256, 0.7485, 0.7538, 0.7584];
err_phonon=[0.0408, 0.0594, 0.0636, 0.0371, 0.0559, 0.0775, 0.078, 0.0623, 0.0767, 0.0629, 0.0218, 0.0474, 0.0219, 0.0385, 0.028, 0.0418, 0.05066, 0.0532, 0.0388, 0.0577, 0.0669, 0.0383, 0.0293];
phase=[0.1042, 0.1802, 0.462, 0.4873, 0.7645, 0.6274, 0.7728, 0.7606, 0.789, 0.9098, 0.2632, 0.3949, 0.4365, 0.5494, 0.6577, 0.7447, 0.8974, 0.0853, 0.0596, 0.1649, 0.2252, 0.275, 0.3286];
err_phase=[0.0124, 0.0168, 0.0167, 0.0108, 0.019, 0.0226, 0.024, 0.0218, 0.0194, 0.0228, 0.0088, 0.0173, 0.0075, 0.0141, 0.0082, 0.0131, 0.0118, 0.0147, 0.014, 0.0212, 0.0283, 0.0161, 0.0118];
secular=0.7676;%MHz

%%%%%%%%%%Calibration
freq_cal=[50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 110., 115., 120., 125., 130., 135., 140., 145., 150., 155., 165.];
ampl_cal=[6., 5.75, 4.25, 8., 32.5, 25., 14., 45., 37.5, 35., 28., 35., 70., 150., 180., 140., 60., 44., 84., 116., 145., 176., 115.];
phonon_cal=[0.8484, 1.033, 1.0339, 0.66, 1.0089, 1.1545, 0.9274, 0.6472, 1.1268, 1.0245, 0.9158, 0.9313, 1.0384, 0.8842, 0.549, 0.5236, 0.7457, 0.7852, 1.0638, 0.7198, 0.7013, 0.7912, 0.6568];
err_phonon_cal=[0.0585, 0.0776, 0.0446, 0.0618, 0.0612, 0.0748, 0.0656, 0.0343, 0.0933, 0.0526, 0.0367, 0.0362, 0.0492, 0.0252, 0.0236, 0.0283, 0.0341, 0.0345, 0.0498, 0.0517, 0.0482, 0.036, 0.0167];
phase_cal=[0.3077, 0.3777, 0.4972, 0.4386, 0.6533, 0.4128, 0.4665, 0.3279, 0.2478, 0.3286, 0.4513, 0.4717, 0.3366, 0.02687, 0.2024, 0.0937, 0.1766, 0.289, 0.2286, 0.2136, 0.293, 0.3627, 0.2747];
err_phase_cal=[0.0214, 0.0217, 0.0111, 0.0252, 0.0161, 0.0179, 0.0186, 0.0155, 0.0264, 0.01575, 0.0106, 0.0101, 0.0144, 0.009, 0.0135, 0.0151, 0.014, 0.0138, 0.0148, 0.0226, 0.0216, 0.0134, 0.0081];

%%%Data sheet
data_freq=[10, 20, 25, 30, 35, 40, 50, 55, 60, 65, 70, 80, 85, 90, 95, 100, 105, 110, 120, 125, 130, 140, 145, 150, 160, 170, 200];
data_trans=[0.15, 0.2, 0.24, 0.24, 0.29, 0.28, 0.27, 0.3, 0.32, 0.28, 0.33, 0.44, 0.48, 0.53, 0.59, 0.76, 1.43, 3.33, 10.34, 14.11, 17.69, 24.23, 27.16, 29.99, 35.25, 39.97, 52.67];
data_time_delay=[7.94, 8, 8.05, 8.18, 8.3, 8.48, 8.78, 9.03, 9.21, 9.5, 9.88, 10.77, 11.44, 12.32, 13.82, 16.17, 18.4, 18.11, 11.27, 8.59, 6.76, 4.69, 4.07, 3.51, 2.84, 2.33, 1.36]*10^-3;%us
data_phase_delay=data_time_delay.*data_freq;
data_phase=[0.0794, 0.16, 0.2013, 0.2454, 0.2905, 0.3392, 0.439, 0.4966, 0.0526, 0.1175, 0.1916, 0.3616, 0.4724, 0.1088, 0.3129, 0.117, 0.432, 0.4921, 0.3524, 0.0737, 0.3788, 0.1566, 0.0902, 0.0265, 0.4544, 0.3961, 0.272];
span=linspace(50,165,1000);
interp_trans=interp1(data_freq,data_trans,span);
interp_phase=interp1(data_freq,data_phase_delay,span);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Constant organization
%e^4/(8hbar m^3 r^6)
const=9.6287968*10^10;%MHz^7
time=1000;%us
kapa=0.714*(1/42);
phi_plus=0.;%turns


coeff=const*cos(2*pi*phi_plus/2)^2*time^2*kapa^2;%MHz^5
scale=coeff./((2*pi*secular).*((2*pi*freq_dipole).^2-(2*pi*secular).^2).^2);

V_sb=.01;
g_cal=(phonon_cal./(scale.*ampl_cal.^2*V_sb^2)).^(1/2);
err_g_cal=(1./(scale.*ampl_cal.^2*V_sb^2)).^(1/2)*1/2.*phonon_cal.^(-1/2).*err_phonon_cal;
g=(phonon./(scale.*ampl_c.^2*V_sb^2)).^(1/2);
err_g=(1./(scale.*ampl_c.^2*V_sb^2)).^(1/2)*1/2.*phonon.^(-1/2).*err_phonon;
db_g=-20*log10(g./g_cal);
err_db=20/log(10)*(err_g./g+err_g_cal./g_cal);

figure();
plot(span, interp_trans,'LineWidth',1.);
hold on
errorbar(freq_dipole,db_g,err_db,'MarkerSize', 8,'LineWidth',1.);
ylabel('Insertion Loss (dB)')
xlabel('Frequency (2\pi*MHz)')
xlim([50,165]);
ylim([-2,40]);

delta_phase=phase-phase_cal;

figure();
plot(span,mod(interp_phase-0.5,1));
hold on
adjust=[-1,-1,-1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,-1,-1,-1,-1,-1,0];
plot(freq_dipole,adjust+mod((phase-phase_cal),1));
plot(span,interp_phase-0.5);

samplex=[0,30,60,90,120,150];
sampley=[0.2,0.5,0.8,0.8,0.5,0.2];
sampledx=[0.01,0.01,0.01,0.01,0.01,0.01];
sampledy=[0.1,0.1,0.1,0.1,0.1,0.1];
fitChiSquare(samplex,sampley,a*sin(b*x+c),[1,1,1],sampledx,sampledy)
