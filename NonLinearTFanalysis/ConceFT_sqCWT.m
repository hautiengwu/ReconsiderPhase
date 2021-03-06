function [tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_CWT(t, x, lowfreq, highfreq, alpha, MT, opts, Smooth, Hemi) ;


N = length(t) ;

tfrsq = zeros(N, round((highfreq-lowfreq)./alpha)) ;



	%% ordinary SST
fprintf(['Run ordinary CWT-SST (Smooth = ',num2str(Smooth),', Hemi = ',num2str(0),')\n']) ;

[tfr, tfrtic, tfrsq, tfrsqtic] = sqCWTbase(t, x, lowfreq, highfreq, alpha, opts, Smooth, 0);


%===========================
    %% get the ConceFT	

ConceFT = tfrsq ;


if MT > 1
	fprintf(['CWT-ConceFT total (Smooth = ',num2str(Smooth),', Hemi = ',num2str(Hemi),'): ',num2str(MT),'; now:     ']) ;

	for ii = 1: MT

	    fprintf('\b\b\b\b') ;   tmp = sprintf('%4d',ii) ; fprintf([tmp]) ;

		[~, ~, tfrsqX, tfrsqtic] = sqCWTbase(t, x, lowfreq, highfreq, alpha, opts, Smooth, Hemi);
	    ConceFT = ConceFT + tfrsqX ;
	
	end
	fprintf('\n') ;

	ConceFT = ConceFT ./ (MT+1) ;
end

end

