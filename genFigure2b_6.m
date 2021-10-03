clear ; close all ;

scrsz = get(0,'ScreenSize');
% from 'With T-P/20161117_35_3176823_TP' ;
% export_fig can be found online
addpath('/Users/hautiengwu/Dropbox/code/NonLinearTFanalysis/tool/export_fig') ;
addpath('/Users/hautiengwu/Dropbox/code/NonLinearTFanalysis/tool') ;
addpath('/Users/hautiengwu/Dropbox/code/NonLinearTFanalysis') ;


tmp = dlmread(['signal_AWF_125Hz.txt'],'\t',1,0) ;


flowTime = tmp(130:end,1) ;
flow = tmp(130:end,2) ;


tmp = dlmread(['signal_AWP_125Hz.txt'],'\t',1,0) ;

Pressure = tmp(130:end,2) ;


tmp = dlmread(['signal_CO2_1_62.5Hz.txt'],'\t',1,0) ;

EtCO2time = tmp(66:end,1) ;
EtCO2 = tmp(66:end,2) ;


tmp = dlmread(['signal_Resp_62.5Hz.txt'],'\t',1,0) ;

impedancetime = tmp(66:end,1) ;
impedance = tmp(66:end,2) ;

EtCO2 = interp1(EtCO2time, EtCO2, flowTime, 'pchip', 'extrap') ;
impedance = interp1(impedancetime, impedance, flowTime, 'pchip', 'extrap') ;

Pressure = Pressure ./ (quantile(Pressure, .99) - quantile(Pressure, .01)) ;
flow = flow ./ (quantile(flow, .99) - quantile(flow, .01)) ;
EtCO2 = EtCO2 ./ (quantile(EtCO2, .99) - quantile(EtCO2, .01)) ;
impedance = impedance ./ (quantile(impedance, .99) - quantile(impedance, .01)) ;

Pressure = Pressure(1:1.5e4) ;
flow = flow(1:1.5e4) ;
EtCO2 = EtCO2(1:1.5e4) ;
impedance = impedance(1:1.5e4) ;



figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/2])
plot(flowTime(1:1.5e4), Pressure,'k', 'linewidth',2) ;
hold on ;
plot(flowTime(1:1.5e4), flow-1, 'r', 'linewidth',2) ;
plot(flowTime(1:1.5e4), EtCO2-3, 'color', [.6 .6 .6],'linewidth',2) ;
plot(flowTime(1:1.5e4), impedance-4, 'b', 'linewidth',2) ;
legend('Pressure', 'Flow', 'EtCO2', 'impedance') ;
axis tight
xlabel('Time (sec)') ;
ylabel('Signal (a.u.)') ;
set(gca,'fontsize',22) ;
export_fig('Figure2b','-transparent','-m3') ;









%% check PS

Pressure = Pressure - mean(Pressure) ;
flow = flow - mean(flow) ;
EtCO2 = EtCO2 - mean(EtCO2) ;
impedance = impedance - mean(impedance) ;

xhatP = fft(Pressure)*(flowTime(2) - flowTime(1)) ;
xhatf = fft(flow)*(flowTime(2) - flowTime(1)) ;
xhate = fft(EtCO2)*(flowTime(2) - flowTime(1)) ;
xhati = fft(impedance)*(flowTime(2) - flowTime(1)) ;

xi = [1:length(xhatP)/2] / length(xhatP) / (flowTime(2) - flowTime(1)) ;

figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/2])
subplot(141) ;
plot(xi, 2*log2(abs(xhatP(2:end/2+1))), 'k', 'linewidth', 2)
xlabel('Frequency (Hz)') ; ylabel('Power spectrum (log2)') ;
axis([0 4 -20 10]) ; set(gca,'fontsize', 24)

subplot(142) ;
plot(xi, 2*log2(abs(xhatf(2:end/2+1))), 'k', 'linewidth', 2)
xlabel('Frequency (Hz)') ; ylabel('Power spectrum (log2)') ;
axis([0 4 -20 10]) ; set(gca,'fontsize', 24)

subplot(143) ;
plot(xi, 2*log2(abs(xhate(2:end/2+1))), 'k', 'linewidth', 2)
xlabel('Frequency (Hz)') ; ylabel('Power spectrum (log2)') ;
axis([0 4 -20 10]) ; set(gca,'fontsize', 24)

subplot(144) ;
plot(xi, 2*log2(abs(xhati(2:end/2+1))), 'k', 'linewidth', 2)
xlabel('Frequency (Hz)') ; ylabel('Power spectrum (log2)') ;
axis([0 4 -20 10]) ; set(gca,'fontsize', 24)





%% get phase (BPF+Hilbert)
PP = bandpass(Pressure,[0.1,0.4], 1./(flowTime(2) - flowTime(1))) ;
FF = bandpass(flow,[0.1,0.4], 1./(flowTime(2) - flowTime(1))) ;
EE = bandpass(EtCO2,[0.1,0.4], 1./(flowTime(2) - flowTime(1))) ;
II = bandpass(impedance,[0.1,0.4], 1./(flowTime(2) - flowTime(1))) ;

phiP = phase(hilbert(PP)) ; phiP = phiP(10*125:1.5e4-10*125) ;
phif = phase(hilbert(FF)) ; phif = phif(10*125:1.5e4-10*125) ;
phie = phase(hilbert(EE)) ; phie = phie(10*125:1.5e4-10*125) ;
phii = phase(hilbert(II)) ; phii = phii(10*125:1.5e4-10*125) ;

figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)*2/3])
subplot(211) 
plot(flowTime(10*125:1.5e4-10*125)-10, mod(phiP, 10*pi), 'k', 'linewidth',2) ;
hold on;
plot(flowTime(10*125:1.5e4-10*125)-10, mod(phif, 10*pi), 'r', 'linewidth', 2) ;
plot(flowTime(10*125:1.5e4-10*125)-10, mod(phie, 10*pi), 'color', [.6 .6 .6], 'linewidth', 2) ;
plot(flowTime(10*125:1.5e4-10*125)-10, mod(phii, 10*pi), 'b', 'linewidth', 2) ;
axis tight ;
legend('Pressure', 'Flow', 'EtCO2', 'impedance') ;
xlabel('Time (sec)') ;
ylabel('phase (angle)') ;
set(gca,'fontsize',22) ;

subplot(212)
IFP = diff(phiP)*125/2/pi ;
plot(flowTime(10*125:25:1.5e4-10*125-1)-10, IFP(1:25:end), 'k', 'linewidth',2) ;
hold on;
IFf = diff(phif)*125/2/pi ;
plot(flowTime(10*125:25:1.5e4-10*125-1)-10, IFf(1:25:end), 'r', 'linewidth', 2) ;
IFe = diff(phie)*125/2/pi ;
plot(flowTime(10*125:25:1.5e4-10*125-1)-10, IFe(1:25:end), 'color', [.6 .6 .6], 'linewidth', 2) ;
IFi = diff(phii)*125/2/pi ;
plot(flowTime(10*125:25:1.5e4-10*125-1)-10, IFi(1:25:end), 'b', 'linewidth', 2) ;
axis tight; axis([0 inf 0 .6])
legend('Pressure', 'Flow', 'EtCO2', 'impedance') ;
xlabel('Time (sec)') ;
ylabel('instantaneous frequency (Hz)') ;
set(gca,'fontsize',22) ;
export_fig('Figure6a','-transparent','-m3') ;







%% run SST

[~, ~, tfrsq1, ~, tfrsqtic] = ConceFT_sqSTFT_C(Pressure-mean(Pressure),...
    0, 0.02, 1e-4, 25, 125*5*8+1, 1, 6, 1, 0, 0) ;

[~, ~, tfrsq2, ~, tfrsqtic] = ConceFT_sqSTFT_C(flow-mean(flow),...
    0, 0.02, 1e-4, 25, 125*5*8+1, 1, 6, 1, 0, 0) ;

[~, ~, tfrsq3, ~, tfrsqtic] = ConceFT_sqSTFT_C(EtCO2-mean(EtCO2),...
    0, 0.02, 1e-4, 25, 125*5*8+1, 1, 6, 1, 0, 0) ;

[~, ~, tfrsq4, ~, tfrsqtic] = ConceFT_sqSTFT_C(impedance-mean(impedance),...
    0, 0.02, 1e-4, 25, 125*5*8+1, 1, 6, 1, 0, 0) ;


figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/3])
subplot(141) 
imageSQ(flowTime(1:25:1.5e4), tfrsqtic*125 , abs(tfrsq1), 0.995) ;
xlabel('Time (sec)') ;
ylabel('Frequency (Hz)') ;
set(gca,'fontsize',22) ; axis([-inf inf 0 .8]) ;

subplot(142) 
imageSQ(flowTime(1:25:1.5e4), tfrsqtic*125, abs(tfrsq2), 0.995) ;
xlabel('Time (sec)') ;
ylabel('Frequency (Hz)') ;
set(gca,'fontsize',22) ; axis([-inf inf 0 .8]) ;

subplot(143) 
imageSQ(flowTime(1:25:1.5e4), tfrsqtic*125, abs(tfrsq3), 0.995) ;
xlabel('Time (sec)') ;
ylabel('Frequency (Hz)') ;
set(gca,'fontsize',22) ; axis([-inf inf 0 .8]) ;

subplot(144) 
imageSQ(flowTime(1:25:1.5e4), tfrsqtic*125, abs(tfrsq4), 0.995) ;
xlabel('Time (sec)') ;
ylabel('Frequency (Hz)') ;
set(gca,'fontsize',22) ; axis([-inf inf 0 .8]) ; colormap(1-gray)

export_fig('Figure6b','-transparent','-m3') ;


%% get phase via SST
[c] = CurveExt_M(abs(tfrsq1(1:40,:)'), 0.5) ;

fundP = zeros(size(tfrsq1,2), 1) ;
fundf = zeros(size(tfrsq2,2), 1) ;
funde = zeros(size(tfrsq1,2), 1) ;
fundi = zeros(size(tfrsq1,2), 1) ;

for jj = 1: length(fundP)
    fundP(jj) = sum(tfrsq1(c(jj)-5: c(jj)+5, jj)) ;
    fundf(jj) = sum(tfrsq2(c(jj)-5: c(jj)+5, jj)) ;
    funde(jj) = sum(tfrsq3(c(jj)-5: c(jj)+5, jj)) ;
    fundi(jj) = sum(tfrsq4(c(jj)-5: c(jj)+5, jj)) ;
end

fundP = interp1(flowTime(1:25:1.5e4), fundP, flowTime(1:1.5e4), 'pchip', 'extrap') ; 
fundP = fundP(10*125:1.5e4-10*125) ;
fundf = interp1(flowTime(1:25:1.5e4), fundf, flowTime(1:1.5e4), 'pchip', 'extrap') ; 
fundf = fundf(10*125:1.5e4-10*125) ;
funde = interp1(flowTime(1:25:1.5e4), funde, flowTime(1:1.5e4), 'pchip', 'extrap') ; 
funde = funde(10*125:1.5e4-10*125) ;
fundi = interp1(flowTime(1:25:1.5e4), fundi, flowTime(1:1.5e4), 'pchip', 'extrap') ; 
fundi = fundi(10*125:1.5e4-10*125) ;
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)*2/3])
subplot(211)
plot(flowTime(10*125:1.5e4-10*125)-10, mod(phase(fundP), 10*pi), 'k', 'linewidth',2) ;
hold on;
plot(flowTime(10*125:1.5e4-10*125)-10, mod(phase(fundf), 10*pi), 'r', 'linewidth', 2) ;
plot(flowTime(10*125:1.5e4-10*125)-10, mod(phase(funde), 10*pi), 'color', [.6 .6 .6], 'linewidth', 2) ;
plot(flowTime(10*125:1.5e4-10*125)-10, mod(phase(fundi), 10*pi), 'b', 'linewidth', 2) ;
axis tight
legend('Pressure', 'Flow', 'EtCO2', 'impedance') ;
xlabel('Time (sec)') ;
ylabel('phase (angle)') ;
set(gca,'fontsize',22) ;
subplot(212)
IFP = diff(phase(fundP))*125/2/pi ;
plot(flowTime(10*125:25:1.5e4-10*125-1)-10, IFP(1:25:end), 'k', 'linewidth',2) ;
hold on;
IFf = diff(phase(fundf))*125/2/pi ;
plot(flowTime(10*125:25:1.5e4-10*125-1)-10, IFf(1:25:end), 'r', 'linewidth', 2) ;
IFe = diff(phase(funde))*125/2/pi ;
plot(flowTime(10*125:25:1.5e4-10*125-1)-10, IFe(1:25:end), 'color', [.6 .6 .6], 'linewidth', 2) ;
IFi = diff(phase(fundi))*125/2/pi ;
plot(flowTime(10*125:25:1.5e4-10*125-1)-10, IFi(1:25:end), 'b', 'linewidth', 2) ;
axis tight
legend('Pressure', 'Flow', 'EtCO2', 'impedance') ;
xlabel('Time (sec)') ;
ylabel('instantaneous frequency (Hz)') ;
set(gca,'fontsize',22) ;
export_fig('Figure6c','-transparent','-m3') ;
