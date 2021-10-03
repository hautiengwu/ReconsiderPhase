function [tfr, tfrtic, tfrsq, tfrsqtic] = sqSTFTmedian(x, lowFreq, highFreq, alpha, HOP, WinLen, dim, supp, MT, Smooth);
%


fprintf(['STFT-ConceFT total (Smooth = ',num2str(Smooth),') + Complex sphere + median\n']) ;




%=======================================
	% prepare for the window
	% h = dim * WinLen
[h0, Dh0, ~] = hermf(WinLen, dim, supp) ;
h0 = h0' ; Dh0 = Dh0' ;

rv = zeros(dim, MT) ;
rv(:, 1) = [1 zeros(1, dim-1)]' ;
for ii = 2: MT
    tmp = randn(dim, 1) + sqrt(-1)*randn(dim, 1) ; 
	rv(:, ii) = tmp ./ norm(tmp) ;
end

h = h0 * rv ;
Dh = Dh0 * rv ;



%=======================================
	% prepare for the data
[xrow,xcol] = size(x) ;
t = [1:length(x)] ;
tLen = length(t(1:HOP:length(x))) ;

	% for tfr
N = length([-0.5+alpha:alpha:0.5]) ;

	% for tfrsq
Lidx = round( (N/2)*(lowFreq/0.5) ) + 1 ; 
Hidx = round( (N/2)*(highFreq/0.5) ) ; 
fLen = Hidx - Lidx + 1 ;



%====================================================================
	%% check input signals
if (xcol~=1),
    error('X must have only one column');
elseif highFreq > 0.5
    error('TopFreq must be a value in [0, 0.5]');
elseif (HOP < 1) | (rem(HOP,1)) 
    error('HOP must be an integer value >= 1');
end; 

[hrow,hcol] = size(h); Lh = (hrow-1)/2; 
if rem(hrow,2)==0
    error('H must be a smoothing window with odd length');
end;


%====================================================================
	%% prepare for the output
tfr = zeros(N/2, tLen); 	% for h
tfrtic = linspace(0, 0.5, N/2)' ;
tfrsq = zeros(fLen, tLen); 
tfrsqtic = linspace(lowFreq, highFreq, fLen)' ;


Ex = mean(abs(x).^2);
Threshold = 1.0e-8*Ex;  % originally it was 1e-6*Ex


%====================================================================
	% prepare for the smoothing step
Mid = round(length(tfrsqtic)/2) ;
Delta = 20*(tfrsqtic(2)-tfrsqtic(1)).^2 ;
weight = exp(-(tfrsqtic(Mid-10:Mid+10)-tfrsqtic(Mid)).^2/Delta) ;
weight = weight ./ sum(weight) ;
weightIDX = [Mid-10:Mid+10] - Mid ;




%====================================================================
	% run ConceFT
for tidx = 1:tLen,

    ti = t((tidx-1)*HOP+1); 
    tau = -min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
    indices= rem(N+tau,N)+1;
	
	for nidx = 1: MT
    	norm_h(nidx) = norm(h(Lh+1+tau, nidx));
	end

	tf0 = zeros(N, MT) ; tf1 = zeros(N, MT) ;
    tf0(indices, :) = repmat(x(ti+tau), 1, MT) .* conj( h(Lh+1+tau, :)) * diag(1./norm_h) ;
    tf1(indices, :) = repmat(x(ti+tau), 1, MT) .* conj(Dh(Lh+1+tau, :)) * diag(1./norm_h) ;
    tf0 = fft(tf0, [], 1) ; tf0 = tf0(1:N/2, :) ;
    tf1 = fft(tf1, [], 1) ; tf1 = tf1(1:N/2, :) ;

		% get the first order omega
	omega = zeros(size(tf1)) ;
	avoid_warn = find(tf0~=0);
	omega(avoid_warn) = imag(N*tf1(avoid_warn)./tf0(avoid_warn)/(2.0*pi)) ;
	omega = round(median(omega, 2)) ;


	sst = zeros(fLen,1) ;

    for jcol = 1: N/2,
  		if abs(tfr(jcol)) > Threshold,

   	    	jcolhat = jcol - omega(jcol) ;

   	    	if (jcolhat <= Hidx) & (jcolhat >= Lidx)

				if Smooth 
					IDXb = find((jcolhat-Lidx+1+weightIDX <= Hidx) & (jcolhat-Lidx+1+weightIDX >= Lidx)) ;
					IDXa = jcolhat-Lidx+1+weightIDX(IDXb) ; 


                   	sst(IDXa) = sst(IDXa) + sum(tf0(jcol, :))*weight(IDXb) ;

				else

	   				sst(jcolhat-Lidx+1) = sst(jcolhat-Lidx+1) + sum(tf0(jcol, :)) ;

				end
        	end



  		end;
    end;


	tfr(:, tidx) = tf0(1:N/2, 1) ;
	tfrsq(:, tidx) = sst ;

end;

