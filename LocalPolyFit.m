function [y,varargout] = LocalPolyFit(xi,yi,x,p,h,varargin)
% LOCALPOLYFIT - Local polynomial regression estimator, or polynomial fit (see Fan and Gijbels section 2.3)
%  Y = LOCALPOLYFIT(XI,YI,X,P,H) implements the local polynomial estimator M(X) = SUM_{N=0}^{N=P} BETA_N*(X-XI)^N 
% by assigning weights within a window of width 2h centered on the point of estimation X. The dimension of Y is M by (P+1) 
% where M is the length of XI and the K-th column is the BETA_(K-1) coefficient of the polynomial
% XI and YI are vectors of the same size of observations of the independent variable X and dependent variable Y. 
% XI is a vector of points where to estimate the function m(X). The Kernel used is the Epanechnikov symmetric beta function
% H is the bandwidth, i.e. the half-width of the weighting window around each point XI where the function is estimated.
% [Y,YE] = LOCALPOLYFIT(XI,YI,X,P,H,E) return the error YE when the rms error for the observations E are specified. 
% E is either a scalar, the same value for all points YI, or is a vector of the same size than YI. 
%The dimension of YE is M by (P+1) where M is the length of X and the K-th column is the error for the BETA_K coefficient of the polynomial

% Shane Elipot, 2014, version 0.2
    
% sort the input?
[xi,I] = sort(xi);
yi = yi(I);

xi = xi(:);
yi = yi(:);

[n1,n2] = size(xi);

% how many optional output arguments?
nout = nargout - 1;

if ~isempty(varargin)
    if length(varargin{:}) == 1
        e = varargin{:}*ones(size(xi));
    elseif (length(varargin{:}) == n1)
        e = varargin{:};
    end
end

% initialize the output
y = NaN*ones(length(x),p+1);
%if nargout==1 & ~isempty(varargin)
%    disp('You provided the input error but did not specify the error output as an output argument?');
%elseif nargout>1 & ~isempty(varargin)
if nout == 1
    ye = NaN*ones(length(x),p+1);
end
    
for k = 1:length(x)

    w = kernelH(xi-x(k),h);

    % need to make sure the problem is not underdetermined: if it is, no estimation possible, or lower the order?
    q = find(w~=0);
    w = w(q);
    w = w/sum(w);

    xi2 = xi(q);
    yi2 = yi(q);

    if length(q)>=p+1
        X = zeros(length(xi2),p+1);
        z = xi2-x(k);
        z = z(:);
        for j = 0:p
            X(:,j+1) = z.^j; 
        end
        W = diag(w,0);
        %W = diag(ones(size(w))/length(w),0);
        Y = yi2(:);        
        A = inv(transpose(X)*W*X)*transpose(X)*W;
        beta = A*Y;
        %C(k) = cond(A);
        y(k,:) = beta.';
     
        if nout == 1
            S = diag(w.^2.*e(q).^2);
            betaV = inv(transpose(X)*W*X)*(transpose(X)*S*X)*inv(transpose(X)*W*X);
            ye(k,:) = [diag(betaV).'].^0.5;
        end
        
    end

end

 if nout == 1
     varargout(1) = {ye};
 end
 
function W = kernelH(t,h);

    g = 1;
    W = kernelB(t/h,g)/h;

function K = kernelB(t,g); % beta family kernel for g = 0,1,2, ...
    
% the kernel is zero outside of the normalized bandwith 1 by construction
    K = (1-t.^2);
    qn = find(K<0);
    K(qn) = 0;
    K = 0.75*K.^g;
    

