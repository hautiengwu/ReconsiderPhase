scrsz = get(0,'ScreenSize');

% export_fig can be found online

initstate(1) ;
t = [1:2048]/128 ; 

x1 = exp(i*2*pi*t) ;

figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
subplot(3,1,1) ;
plot(t, real(x1), 'k', 'linewidth', 2)
hold on ;
plot(t, imag(x1), 'r', 'linewidth', 2)
xlabel('Time (Sec)') ; ylabel('a.u.') ;
set(gca, 'fontsize', 24) ;

A = 1+0.3*cos(t).*(smooth(rand(size(t))-1, 256, 'loess')') ;
x2 = A.*exp(i*2*pi*(t+0.02*t.^2.5)) ;
subplot(3,1,2) ;
plot(t, real(x2), 'k', 'linewidth', 2)
hold on ;
plot(t, imag(x2), 'r', 'linewidth', 2)
xlabel('Time (Sec)') ; ylabel('a.u.') ;
set(gca, 'fontsize', 24) ;
 
x3 = exp(i*2*pi*(t+0.02*t.^2.5)) + 0.95*exp(i*(4*pi*(t+0.02*t.^2.5+3/4/pi))) ;
subplot(3,1,3) ;
plot(t, real(x3), 'k', 'linewidth', 2)
hold on ;
plot(t, imag(x3), 'r', 'linewidth', 2)
xlabel('Time (Sec)') ; ylabel('a.u.') ;
set(gca, 'fontsize', 24) ;

export_fig('Figure1','-transparent','-pdf') ;




%% 
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])

tt = [1-256:2048+256]/128 ; 

x3 = exp(i*2*pi*(tt+0.02*tt.^2.5)) + 0.95*exp(i*(4*pi*(tt+0.02*tt.^2.5+3/4/pi))) ;
x4 = exp(i*2*pi*(tt+0.02*tt.^2.5)) + 1.01*exp(i*(4*pi*(tt+0.02*tt.^2.5+3/4/pi))) ;


[~, ~, tfrsq1, ~, tfrsqtic] = ConceFT_sqSTFT_C(real(x3)'-mean(real(x3)),...
    0, 0.1, 1e-4, 1, 128*6+1, 1, 6, 1, 1, 0) ;

[~, ~, tfrsq2, ~, tfrsqtic] = ConceFT_sqSTFT_C(real(x4)'-mean(real(x4)),...
    0, 0.1, 1e-4, 1, 128*6+1, 1, 6, 1, 1, 0) ;

tfrsq1 = tfrsq1(:, 256+1: end-256) ;
tfrsq2 = tfrsq2(:, 256+1: end-256) ;
x3 = x3(256+1: end-256) ;
x4 = x4(256+1: end-256) ;

c = ceil((1+0.02*2.5.*t.^1.5)./ (128*(tfrsqtic(2)-tfrsqtic(1)))) ;

recon1 = [] ; recon2 = [] ;
for jj = 1:2048
   recon1(jj) = sum(tfrsq1(c(jj)- 8: c(jj)+8, jj)) ;
   recon2(jj) = sum(tfrsq2(c(jj)- 8: c(jj)+8, jj)) ;
end


subplot(241)
plot(t(1:1024), x3(1:1024), 'k','linewidth', 2) ; axis tight ;
xlabel('Time (Sec)') ; ylabel('a.u.') ;
set(gca, 'fontsize', 24) ;

subplot(242)
plot(t, 2*pi*(t+0.02*t.^2.5)-phase(hilbert(real(x3))), 'r', 'linewidth', 2) ; hold on ;
%plot(t, phase(hilbert(real(x3))), 'b--', 'linewidth', 2) ;
xlabel('Time (Sec)') ; ylabel('phase diffence (angle)') ;
set(gca, 'fontsize', 24) ; axis([-inf inf -pi pi]) ; %axis([-inf inf 0 350]) ;

subplot(243)
imageSQ(t, tfrsqtic*128 , abs(tfrsq1), 0.995) ;
xlabel('Time (sec)') ;
ylabel('Frequency (Hz)') ;
set(gca,'fontsize',24) ; axis([-inf inf 0 5]) ;

subplot(244)
plot(t, 2*pi*(t+0.02*t.^2.5)-phase(recon1), 'b', 'linewidth', 2) ; hold on ;
%plot(t, phase(recon1), 'k--', 'linewidth', 2) ;
xlabel('Time (Sec)') ; ylabel('phase diffence (angle)') ; axis([-inf inf -pi pi]) ;
set(gca, 'fontsize', 24) ;  


subplot(245)
plot(t(1:1024), x4(1:1024), 'k','linewidth', 2) ; axis tight ;
xlabel('Time (Sec)') ; ylabel('a.u.') ;
set(gca, 'fontsize', 24) ;

subplot(246)
plot(t, 2*pi*(t+0.02*t.^2.5)-phase(hilbert(real(x4))), 'r', 'linewidth', 2) ; hold on ;
%plot(t, phase(hilbert(real(x4))), 'b--', 'linewidth', 2) ;
xlabel('Time (Sec)') ; ylabel('phase diffence (angle)') ;
set(gca, 'fontsize', 24) ; %axis([-inf inf 0 350]) ;


subplot(247)
imageSQ(t, tfrsqtic*128 , abs(tfrsq2), 0.995) ;
xlabel('Time (sec)') ;
ylabel('Frequency (Hz)') ;
set(gca,'fontsize',24) ; axis([-inf inf 0 5]) ;

subplot(248)
plot(t, 2*pi*(t+0.02*t.^2.5)-phase(recon2), 'b', 'linewidth', 2) ; hold on ;
%plot(t, phase(recon2), 'k--', 'linewidth', 2) ;
xlabel('Time (Sec)') ; ylabel('phase diffence (angle)') ; axis([-inf inf -pi pi]) ;
set(gca, 'fontsize', 24) ;

export_fig('Figure2','-transparent','-pdf') ;
