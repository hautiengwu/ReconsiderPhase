function imageCWT(t, ytic, image, type, yN)
%
% Display time-frequency result of the Synchrosqueezing transform ver 0.2
% ImageSQ(t, ytic, image, 'SQ', yN)
% by Hau-tieng Wu 2011-06-20 (hauwu@math.princeton.edu)
% by Hau-tieng Wu 2012-03-03 (hauwu@math.princeton.edu)
% by Hau-tieng Wu 2012-12-23 (hauwu@math.princeton.edu)
%

if nargin < 5 ; yN = 5 ; end


fz = 18;

[n, nscale] = size(image);
xtic = t;

QQ = quantile(image(:), 0.999);
image(find(image>QQ)) = QQ; clear QQ;

imagesc(xtic, 1:nscale, image'); axis xy;
set(gca,'YTick',floor(nscale/yN):floor(nscale/yN):nscale);

if strcmp(type, 'SQ')
    set(gca,'YTickLabel',round(100*ytic(get(gca,'YTick')))/100);
elseif strcmp(type, 'CWT')
    set(gca,'YTickLabel',round(100*(t(end)-t(1))./ytic(get(gca,'ytick')))/100);
end

set(gca,'fontsize',fz);
colormap(1-gray);

end




