function yhat = sgolay(y,order, halfWidth)
% Savitzky-Golay filter
% inputs:
%   y - data to be smoothed (N x L)
%   order - order of local polynomial to fit (suggest 3)
%   halfWidth - range of data to smooth over

L = size(y,2);
yhat = zeros(size(y));

% check for cached coefficients to speed up computation
%{
H_fstr = ['H_',num2str(order),'_',num2str(halfWidth),'.mat'];
if(exist(H_fstr))
    eval(['load ',H_fstr]);
    disp('loaded')
else
  
%}
% compute coefficients
x = [-halfWidth:halfWidth];
X = ones(size(x));
for k=1:order
    X = [ones(size(x)); X.*repmat(x,size(X,1),1)];
end

Bhat = pinv(X);
H = Bhat*X; % symmetric matrix containing coefficients
hh = H(:,halfWidth+1);
%eval(['save ',H_fstr,' H']);
%disp('saved')
%end
for j=1:size(y,1)
    %for i=halfWidth+1:L-halfWidth % middle section
    %    yhat(j,i) = y(j,i-halfWidth:i+halfWidth)*hh; % replace this line with a 'conv'
    %end
    yhat(j,:) = conv(y(j,:),hh,'same');
    
    % smooth at edges with a full width
    for i=1:halfWidth
        yhat(j,i) = y(j,1:2*halfWidth+1)*H(:,i);
        yhat(j,end-i+1) = y(j,end-2*halfWidth:end)*H(:,end-i+1);
    end
end
