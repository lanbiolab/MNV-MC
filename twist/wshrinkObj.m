function [x,objV] = wshrinkObj(x,rho,sX, isWeight,mode)

if isWeight == 1
    C = sqrt(sX(3)*sX(2));
end
if ~exist('mode','var')
    % mode = 1是采用lateral slice的方法
    % mode = 2是采用front slice的方法
    % mode = 3是采用top slice的方法
    mode = 1;
end

X=reshape(x,sX);%传进来的是一长串，故先变回sX里保存的三维
if mode == 1
    Y=X2Yi(X,3);
elseif mode == 3
    Y=shiftdim(X, 1);%维度左移一位NNK变成NKN
else
    Y = X;
end
Yhat = fft(Y,[],3);%傅里叶变换

objV = 0;
if mode == 1
    n3 = sX(2);
elseif mode == 3
    n3 = sX(1);
else%这里决定算法1里的n3是谁，mode3时是NKN，故n3是N
    n3 = sX(3);%mode2时是NNK,故n3是K
end

endValue = int16(n3/2+1);%n3只需要算一半
for i = 1:endValue
    [uhat,shat,vhat] = svd(full(Yhat(:,:,i)),'econ');
    if isWeight
        weight = C./(diag(shat) + eps);
        tau = rho*weight;
        shat = soft(shat,diag(tau));
    else
        tau = rho;
        shat = max(shat - tau,0);
    end
    
    objV = objV + sum(shat(:));
    Yhat(:,:,i) = uhat*shat*vhat';
    if i > 1
        Yhat(:,:,n3-i+2) = conj(uhat)*shat*conj(vhat)';
        objV = objV + sum(shat(:));
    end
end

if isinteger(n3/2)
    [uhat,shat,vhat] = svd(full(Yhat(:,:,endValue+1)),'econ');
    if isWeight
        weight = C./(diag(shat) + eps);
        tau = rho*weight;
        shat = soft(shat,diag(tau));
    else
        tau = rho;
        shat = max(shat - tau,0);
    end
    objV = objV + sum(shat(:));
    Yhat(:,:,endValue+1) = uhat*shat*vhat';
end

Y = ifft(Yhat,[],3);
if mode == 1
    X = Yi2X(Y,3);
elseif mode == 3
    X = shiftdim(Y, 2);
else
    X = Y;
end

x = X(:);

end
