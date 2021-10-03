clear ; close all ;

scrsz = get(0,'ScreenSize');
folder1 = ['/Users/hautiengwu/Dropbox/___working tex 20191219/'...
    '_WORKING_KirkAymen LBNP/LBNP code/database/high tolerance'] ;

folder2 = ['/Users/hautiengwu/Dropbox/___working tex 20191219/'...
    '_WORKING_KirkAymen LBNP/LBNP code/database/low tolerance'] ;


addpath('/Users/hautiengwu/Dropbox/code/NonLinearTFanalysis/tool/export_fig') ;
addpath('/Users/hautiengwu/Dropbox/code/NonLinearTFanalysis/tool') ;
addpath('/Users/hautiengwu/Dropbox/code/NonLinearTFanalysis') ;


ggg1=dir([folder1,'/*.mat']) ;

for jj = length(ggg1):-1:1
    if length(strfind(ggg1(jj).name, 'EKG')) == 0
        ggg1(jj) = [] ;
    end
end



ggg2 = dir([folder2,'/*.mat']) ;

for jj = length(ggg2):-1:1
    if length(strfind(ggg2(jj).name, 'EKG')) == 0
        ggg2(jj) = [] ;
    end
end




KKK = 1 ;

HOP = 10 ;

for caseNO = 1: length(ggg1) + length(ggg2)
    
    
    if caseNO <= length(ggg1)
        load([folder1,'/', ggg1(caseNO).name]) ;
    else
        load([folder2,'/', ggg2(caseNO-length(ggg1)).name]) ;
    end
    
    
    
    ECG = data(datastart(1): dataend(1)) ;
    PVP0 = data(datastart(2): dataend(2)) ;
    PPG = data(datastart(3): dataend(3)) ;
    
    % at least 5 mins
    if length(ECG)/100 > 5*60
        
        
        R = HRCFTG(ECG, 100) ;
        IHRt = interp1(R(2:end)/100, 100./diff(R), [1: 100/4: length(ECG)]/100, 'pchip', 'extrap') ;
        IHRt = medfilt1(IHRt, 4) ;

        
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
        
        
        
        
        
        
        
        
        %% get phase (BPF+Hilbert)
        PVPb = bandpass(PVP,[0.5,2], 100) ;
        PPGb = bandpass(PPG,[0.5,2], 100) ;
        
        phiPPG = phase(hilbert(PVPb)) ;
        phiPVP = phase(hilbert(PPGb)) ;
        
        
        IHRppg0 = diff(phiPPG)*10/2/pi ; IHRppg0=[IHRppg0 IHRppg0(end)]; 
        IHRpvp0 = diff(phiPVP)*10/2/pi ; IHRpvp0=[IHRpvp0 IHRpvp0(end)]; 
        
        IHRppg0 = medfilt1(IHRppg0, 100) ;
        IHRpvp0 = medfilt1(IHRpvp0, 100) ;
        
        IHRppg = IHRppg0(1:25:end) ;
        IHRpvp = IHRpvp0(1:25:end) ;
        
        IHRallMedian_PPGECG_H(KKK) = median(IHRt - IHRppg) ;
        IHRallMedian_PVPECG_H(KKK) = median(IHRt - IHRpvp) ;
        IHRall_PPG_H{KKK} = IHRppg ;
        IHRall_PVP_H{KKK} = IHRpvp ;
        
        
        
        
        %% run SST
        [~, ~, tfrsq1, ~, tfrsqtic] = ConceFT_sqSTFT_C(PPG'-mean(PPG),...
            0, 0.1, 1e-4, HOP, 100*8+1, 1, 6, 1, 0, 0) ;
        
        [~, ~, tfrsq2, ~, tfrsqtic] = ConceFT_sqSTFT_C(PVP'-mean(PVP),...
            0, 0.1, 1e-4, HOP, 100*8+1, 1, 6, 1, 0, 0) ;
        
        
        
        
        
        %% get phase via SST
        [c] = CurveExt_M(abs(tfrsq1(1:300,:)'), 0.5) ;
        
        fundPPG = zeros(size(tfrsq1,2), 1) ;
        fundPVP = zeros(size(tfrsq2,2), 1) ;
        
        for jj = 1: length(fundPPG)
            fundPPG(jj) = sum(tfrsq1(max(1, c(jj)-10): c(jj)+10, jj)) ;
            fundPVP(jj) = sum(tfrsq2(max(1, c(jj)-10): c(jj)+10, jj)) ;
        end
        
        fundPPG = fundPPG(10: end-10) ;
        fundPVP = fundPVP(10: end-10) ;
        
        
        IHRppg0 = diff(phase(fundPPG))*10/2/pi ;
        IHRpvp0 = diff(phase(fundPVP))*10/2/pi ;
        IHRppg0 = medfilt1(IHRppg0, 10) ;
        IHRpvp0 = medfilt1(IHRpvp0, 10) ;
        
        IHRppg = interp1(time(10*HOP:HOP:end-10*HOP-1), IHRppg0, ...
            time(1:25:end), 'pchip', 'extrap');
        IHRpvp = interp1(time(10*HOP:HOP:end-10*HOP-1), IHRpvp0, ...
            time(1:25:end), 'pchip', 'extrap');
        
        
        IHRallMedian_PPGECG_S(KKK) = median(IHRt - IHRppg) ;
        IHRallMedian_PVPECG_S(KKK) = median(IHRt - IHRpvp) ;
        IHRall_PPG_S{KKK} = IHRppg ;
        IHRall_PVP_S{KKK} = IHRpvp ;
        IHRall_ECG{KKK} = IHRt ;
        
        KKK = KKK + 1 ;

    end
end


%%
NN = 0 ;
for jj = 1: length(IHRall_PVP_S)
    NN = NN + length(IHRall_PVP_S{jj}) ;
end

PPGallS = zeros(NN, 1) ;
PVPallS = zeros(NN, 1) ;
PPGallH = zeros(NN, 1) ;
PVPallH = zeros(NN, 1) ;
ECGall = zeros(NN, 1) ;


piv = 1 ;
for jj = 1: length(IHRall_PVP_S)
    PPGallS(piv: piv + length(IHRall_PVP_S{jj}) - 1) = IHRall_PPG_S{jj} ;
    PVPallS(piv: piv + length(IHRall_PVP_S{jj}) - 1) = IHRall_PVP_S{jj} ;
    PPGallH(piv: piv + length(IHRall_PVP_S{jj}) - 1) = IHRall_PPG_H{jj} ;
    PVPallH(piv: piv + length(IHRall_PVP_S{jj}) - 1) = IHRall_PVP_H{jj} ;
    ECGall(piv: piv + length(IHRall_ECG{jj}) - 1) = IHRall_ECG{jj} ;
    piv = piv + length(IHRall_PVP_S{jj}) ;
end

save QQ PPGallS PVPallS PPGallH PVPallH