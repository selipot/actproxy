function x0=findzero(x,y)
% FINDZERO - find a zero crossing in a vector
% find zero(s) of y(x)
% employs linear interpolation 
x=x(:);
y=y(:);

%eliminate the gaps
id=find(~isnan(y) & ~isnan(x));
x=x(id);
y=y(id);

if isempty(x) x0=[];return;end;
a=diff(y)./diff(x);
b=-a.*x(1:end-1)+y(1:end-1);

a(a==0)=nan;

xn=-b./a;
inc_mask=(diff(x)>=0);
dec_mask=(diff(x)<0);

id=find((inc_mask & xn>=x(1:end-1) & xn<x(2:end)) | (dec_mask & xn<=x(1:end-1) & xn>x(2:end)));

x0=xn(id);
