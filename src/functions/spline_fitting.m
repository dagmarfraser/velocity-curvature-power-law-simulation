function y_=spline_fitting(x,x_,y,NN,mean_off)

if nargin==5 && mean_off
    y = y - repmat(mean(y,1),[size(y,1),1]);
end

y_pp = splinefit(x,y',NN);
y_   = ppval(y_pp, x_)';
end