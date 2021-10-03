clear ; close all ;

scrsz = get(0,'ScreenSize');
% from high tolerance
% export_fig can be found online

addpath('/Users/hautiengwu/Dropbox/code/NonLinearTFanalysis/tool/export_fig') ;
addpath('/Users/hautiengwu/Dropbox/code/NonLinearTFanalysis/tool') ;
addpath('/Users/hautiengwu/Dropbox/code/NonLinearTFanalysis') ;

load(['Case #31_FINGER Nellcor_EKG.mat']) ;



HOP = 10 ;

ECG = data(datastart(1): dataend(1)) ;
PVP0 = data(datastart(2): dataend(2)) ;
PPG = data(datastart(3): dataend(3)) ;


R = HRCFTG(ECG, 100) ;
IHRt = interp1(R(2:end)/100, 100./diff(R), [1: 100/4: length(ECG)]/100, 'pchip', 'extrap') ;

trend = [] ;
for jj = 1:length(PVP0)
    idx = max(1, jj-100*2): min(length(PVP0), jj+100*2) ;
    trend(jj) = median(PVP0(idx)) ;
end
trend = smooth(trend, 400, 'loess') ;
PVP = PVP0 - trend' ;


PPG = PPG ./ (quantile(PPG, .99) - quantile(PPG, .01)) ;
PVP = PVP ./ (quantile(PVP, .99) - quantile(PVP, .01)) ;

PPG = PPG - mean(PPG) ;
PVP = PVP - mean(PVP) ;
time = [1:length(PPG)] / 100 ;


figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)/2])
plot(time(1:5e3), PPG(com(3,3)+5e3+1:com(3,3)+1e4),'k', 'linewidth',2) ;
hold on ;
plot(time(1:5e3), PPG(com(8,3)+5e3+1:com(8,3)+1e4)-1, 'b', 'linewidth',2) ;
plot(time(1:5e3), PVP(com(3,3)+5e3+1:com(3,3)+1e4)-2, 'r', 'linewidth',2) ;
plot(time(1:5e3), PVP(com(8,3)+5e3+1:com(8,3)+1e4)-3, 'm', 'linewidth',2) ;
legend('PPG (baseline)', 'PPG (-75mmHG)', 'PVP (baseline)', 'PVP (-75mmHG)') ;
axis tight
xlabel('Time (sec)') ;
ylabel('Signal (a.u.)') ;
set(gca,'fontsize',22) ;
export_fig('Figure2c','-transparent','-m2') ;




%% check PS

xhatPPG = fft(PPG)*100 ;
xhatPVP = fft(PVP)*100 ;

xi = 100*[1:length(xhatPPG)/2] / length(xhatPPG) ;

figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)/3])
subplot(121) ;
plot(xi, 2*log2(abs(xhatPPG(2:end/2+1))), 'k', 'linewidth', 2)
xlabel('Frequency (Hz)') ; ylabel('Power spectrum (log2)') ;
axis([0 4 -5 50]) ; set(gca,'fontsize', 24)

subplot(122) ;
plot(xi, 2*log2(abs(xhatPVP(2:end/2+1))), 'k', 'linewidth', 2)
xlabel('Frequency (Hz)') ; ylabel('Power spectrum (log2)') ;
axis([0 4 -5 50]) ; set(gca,'fontsize', 24)






%% get phase (BPF+Hilbert)
PVPb = bandpass(PVP,[0.5,2], 100) ;
PPGb = bandpass(PPG,[0.5,2], 100) ;

phiPPG = phase(hilbert(PVPb)) ;
phiPVP = phase(hilbert(PPGb)) ;

figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)*2/3])
subplot(211) 
plot(time(1:end-1), diff(phiPPG)/(1/100)/2/pi, 'k', 'linewidth',2) ;
hold on;
plot(time(1:end-1), diff(phiPVP)/(1/100)/2/pi, 'r', 'linewidth', 2) ;
axis tight
legend('PPG', 'PVP') ;
xlabel('Time (sec)') ;
ylabel('instantaneous frequency (Hz)') ;
set(gca,'fontsize',22) ; axis([0 inf 0 5])



PVPb = bandpass(PVP,[0.5,4], 100) ;
PPGb = bandpass(PPG,[0.5,4], 100) ;

phiPPG = phase(hilbert(PVPb)) ;
phiPVP = phase(hilbert(PPGb)) ;

subplot(212)
plot(time(1:end-1), diff(phiPPG)/(1/100)/2/pi, 'k', 'linewidth',2) ;
hold on;
plot(time(1:end-1), diff(phiPVP)/(1/100)/2/pi, 'r', 'linewidth', 2) ;
axis tight
legend('PPG', 'PVP') ;
xlabel('Time (sec)') ;
ylabel('instantaneous frequency (Hz)') ;
set(gca,'fontsize',22) ; axis([0 inf 0 5])






export_fig('Figure7a','-transparent','-pdf') ;



%% run SST
[~, ~, tfrsq1, ~, tfrsqtic] = ConceFT_sqSTFT_C(PPG'-mean(PPG),...
    0, 0.1, 1e-4, HOP, 100*8+1, 1, 6, 1, 0, 0) ;

[~, ~, tfrsq2, ~, tfrsqtic] = ConceFT_sqSTFT_C(PVP'-mean(PVP),...
    0, 0.1, 1e-4, HOP, 100*8+1, 1, 6, 1, 0, 0) ;



figure('Position',[1 scrsz(4) scrsz(3)/2 scrsz(4)/3])
subplot(121) 
imageSQ(time(1:HOP:end), tfrsqtic*100 , abs(tfrsq1), 0.995) ;
xlabel('Time (sec)') ;
ylabel('Frequency (Hz)') ;
set(gca,'fontsize',22) ; axis([-inf inf 0 5]) ;

subplot(122) 
imageSQ(time(1:HOP:end), tfrsqtic*100, log(1+abs(tfrsq2)), 0.995) ;
xlabel('Time (sec)') ;
ylabel('Frequency (Hz)') ;
set(gca,'fontsize',22) ; axis([-inf inf 0 5]) ;


export_fig('Figure7b','-transparent','-m2') ;


%% get phase via SST
[c] = CurveExt_M(abs(tfrsq1(1:300,:)'), 0.5) ;

fundPPG = zeros(size(tfrsq1,2), 1) ;
fundPVP = zeros(size(tfrsq2,2), 1) ;

for jj = 1: length(fundPPG)
    fundPPG(jj) = sum(tfrsq1(c(jj)-10: c(jj)+10, jj)) ;
    fundPVP(jj) = sum(tfrsq2(c(jj)-10: c(jj)+10, jj)) ;
end

fundPPG = fundPPG(10: end-10) ;
fundPVP = fundPVP(10: end-10) ;
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)*2/3])
subplot(211)
plot(time(10*HOP:HOP:end-10*HOP), phase(fundPVP), 'r', 'linewidth', 2) ;
hold on;
plot(time(10*HOP:HOP:end-10*HOP), phase(fundPPG), 'k', 'linewidth',2) ;
axis tight
legend('PVP', 'PPG') ;
xlabel('Time (sec)') ;
ylabel('phase (angle)') ;
set(gca,'fontsize',22) ;
subplot(212)
plot(time(10*HOP:HOP:end-10*HOP-1), diff(phase(fundPVP))*10/2/pi, 'r', 'linewidth', 2) ;
hold on;
plot(time(10*HOP:HOP:end-10*HOP-1), diff(phase(fundPPG))*10/2/pi, 'k', 'linewidth',2) ;
axis tight
legend('PVP', 'PPG') ;
xlabel('Time (sec)') ;
ylabel('instantaneous frequency (Hz)') ;
set(gca,'fontsize',22) ;  axis([0 inf 0 5])

export_fig('Figure7c','-transparent','-m3') ;


IHRppg = interp1(time(10*HOP:HOP:end-10*HOP-1), diff(phase(fundPPG))*10/2/pi, ...
    time(1:25:end), 'pchip', 'extrap'); 
IHRpvp = interp1(time(10*HOP:HOP:end-10*HOP-1), diff(phase(fundPVP))*10/2/pi, ...
    time(1:25:end), 'pchip', 'extrap'); 


IHRwholeMedian_PPGECG(caseNO) = median(IHRt - IHRppg) ;
IHRwholeMedian_PVPECG(caseNO) = median(IHRt - IHRpvp) ;
IHRwholeMAD_PPGECG(caseNO) = mad(IHRt - IHRppg) ;
IHRwholeMAD_PVPECG(caseNO) = mad(IHRt - IHRpvp) ;
