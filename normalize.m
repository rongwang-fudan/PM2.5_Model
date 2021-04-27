function xm = normalize(x)
% divided by the std and reduced by the mean
i2=size(x,2);
xm=x;
for i=1:i2
%     xm(:,i) = xm(:,i)./std(xm(:,i),1);
    xm(:,i) = xm(:,i)-mean(xm(:,i),1);
end
return

