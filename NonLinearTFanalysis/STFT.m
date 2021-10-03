function [tfr, tfrtic] = STFT(x, tDS, alpha, h);
%
%	computes the STFT  (modified from TFtoolbox by Hau-tieng Wu)
%
%   Example:
%
%	Hz = 32 ; t=[1/Hz:1/Hz:32]' ;
%	x=cos(2*pi*(4*t+cos(t/2))) ;
%	[h, Dh] = hermf(71, 1, 6) ;
%		%% get the TF representation of the signal x with the frequency range
%		%% [0,1, 0.4]*Hz with the frequency resolution 0.001*Hz
%	[tfr, tfrtic] = STFT(x, 1, 0.01,  h');
%	imageSQ(t, tfrsqtic*Hz, abs(tfrsq), .998) ;
%
%	X     : signal.
%	tDS   : the time axis in the final TF representation is downsampled by tDS (set it to 1 unless you know what you are doing)
%	alpha : the frequency axis resolution
%	H     : frequency smoothing window, H(0) being forced to 1
%	TFR   : STFT
%
%	F. Auger, May-July 1994, July 1995.
%	Copyright (c) 1996 by CNRS (France).
%
%	------------------- CONFIDENTIAL PROGRAM -------------------- 
%	This program can not be used without the authorization of its
%	author(s). For any comment or bug report, please send e-mail to 
%	f.auger@ieee.org 

if nargin < 8
	ODD = 0 ;
end




[xrow,xcol] = size(x) ;
t = [1:length(x)] ;
tLen = length(t(1:tDS:length(x))) ;

	% for tfr
N = length([-0.5+alpha:alpha:0.5]) ;


%====================================================================
	%% check input signals
if (xcol~=1),
    error('X must have only one column');
elseif highFreq > 0.5
    error('TopFreq must be a value in [0, 0.5]');
elseif (tDS < 1) | (rem(tDS,1)) 
    error('tDS must be an integer value >= 1');
end; 

[hrow,hcol] = size(h); Lh = (hrow-1)/2; 
if (hcol~=1)|(rem(hrow,2)==0),
    error('H must be a smoothing window with odd length');
end;


%====================================================================
	%% run STFT and reassignment rule
tfr = zeros(N/2, tLen); 	% for h
tfrtic = linspace(0, 0.5, N/2)' ;



for tidx = 1:tLen,

    ti = t((tidx-1)*tDS+1); 
    tau = -min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
    indices= rem(N+tau,N)+1 ;
    norm_h = norm(h(Lh+1+tau)) ;

	tf0 = zeros(N, 1) ; 
    tf0(indices) = x(ti+tau).*conj( h(Lh+1+tau)) /norm_h;
    tf0 = fft(tf0) ; tfr(:, tidx) = tf0(1:N/2) ;

end;

