%% in this simulation example, how to run ConceFT is demonstrated with the conceFT for the synchrosqueezed STFT. The same flow could be appied to conceFT for the synchrosqueezed CWT or others. Please see the code for details.

addpath('./tool') ;

	%% the sampling rate for the simulated signal
Hz = 32 ;
	%% the sampling time of the simulated signal
time = [1/Hz:1/Hz:16]' ;
	%% the number of sampling points of the simulated signal
N = length(time) ;

    %% fix the random seed for the reproducibility issue
initstate(1) ;


	%% the amplitude modulation of the simulated signal
	%% simulate 1 oscillatory components with dynamics
am1 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
am1 = 2 + am1 ./ max(abs(am1)) ;

    %% the instantaneous frequency of the simulated signal
if1 = smooth(cumsum(randn(N,1)) ./ Hz, 400, 'loess') ;
if1 = 4 + 1 * if1 ./ max(abs(if1)) ;
phi1 = cumsum(if1) / Hz ; 

	%% the simulated signal.
s1 = am1 .* cos(2*pi*phi1) ; 
clean = s1 ;


	%% add noise 
sigma = 0 ;
noise = random('T',4,N,1) ;
noise = sigma * noise ; 
var(noise)
snrdb = 20 * log10(std(clean)./std(noise)) ;
fprintf(['snrdb = ',num2str(snrdb),'\n']) ;

	%% simulated observed time series
xm = clean + noise ;

	%% setup parameters for the SST or ConceFT

	%% the window length. Ideally, it should be chosen so that
	%% roughly 7-10 oscillations (ignore the multiples) are 
	%% included in the window.
WindowLength = 377 ;
	%% this is the bandwith of the chosen window. See hermf.m
	%% in the attached code for details.
WindowBandwidth = 10 ;
SamplingRate = Hz ;
    %% Setup the frequency range for the analysis
	%% The allowed range is 0-0.5
	%% This part might be tricky. This 0-0.5 limitation is 
	%% setup under the assumption that the sampling rate is 1Hz
	%% After the analysis, the sampling rate should be adjusted
	%% so that the result is associated with the original 
	%% sampling rate. 
	%% In this example, the true range is [0, 0.5]*SamplingRate
HighFrequencyLimit = 0.5 ;
LowFrequencyLimit = 0 ;
	%% the frequency axis resolution in the final time-frequency representation
FrequencyAxisResolution = 0.001 ;

	%% pad the signal to alleviate the boundary effect
[xpad, nup, n1, n2] = padsignal(xm, 1) ;

	%% synchrosqueezed short time Fourier transform (STFT)
	%% Output:
	%% tfr: STFT result
	%% tfrtic: frequency axis tic for the STFT
	%% tfrsq: synchrosqueezed STFT (it is equivalent to running ConceFT only one time)
	%% tfrsqtic: frequency axis tic for the tfrsq and ConceFT
[h, Dh] = hermf(WindowLength, 1, WindowBandwidth) ;
[tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFTbase(xm, LowFrequencyLimit, HighFrequencyLimit, FrequencyAxisResolution, 1, h', Dh', 0, 0);

	%% resize the TF representation
tfr = tfr(:, n1+1:n1+length(xm)) ;
tfrsq = tfrsq(:, n1+1:n1+length(xm)) ;

	%% plot the time frequency representation determined by
	%% ConceFT. .995 is the quantile truncation set to avoid 
	%% possible outliers in the final analysis result.
	%% see the code for details.
figure ;
subplot(121)
imageSQ(t, tfrsqtic*SamplingRate, abs(tfr), .999) ; colormap(1-gray) ; title('STFT') ;
subplot(122)
imageSQ(t, tfrsqtic*SamplingRate, abs(tfrsq), .999) ; colormap(1-gray) ; title('SST') ;


