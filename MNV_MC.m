clear;
addpath(genpath(pwd))
load('yale_mtv');
%%
lambda=0.08;%%%%%%%%%%%%%%%%%%%%%%%%%%%0.01~0.2
cls_num = length(unique(Y));
K = length(X); N = size(X{1},2);
for v=1:K
    [X{v}]=NormalizeData(X{v});
end
for k=1:K  
    Z{k} = zeros(N,N); 
    Q1{k} = zeros(N,N);
    G{k} = zeros(N,N); 
    J{k} = zeros(N,N);
    Q3{k} = zeros(N,N);
    E{k} = zeros(size(X{k},1),N);
    Q2{k} = zeros(size(X{k},1),N); 
end
sX = [N, N, K];
%set Default
Isconverg = 0;epson = 0.00001;
mu = 10e-5; max_mu = 10e10; pho_mu = 1.4;
iter = 0;
for k=1:K
    tmp_inv{k} = inv(2*eye(N,N)+X{k}'*X{k});
end

while(Isconverg == 0)
    fprintf('----processing iter %d--------\n', iter+1);
    for k=1:K
        %-------------------1 update Z^k-------------------------------
        tmp = (X{k}'*Q2{k} - Q3{k} - Q1{k} )/mu +G{k}+J{k}+ X{k}'*X{k} - X{k}'*E{k};
        Z{k}=tmp_inv{k}*tmp;
    end
    %-------------------2 update E^k-------------------------------
    F = [];
    for k=1:K
        tmp = X{k}-X{k}*Z{k}+Q2{k}/mu;
        F = [F;tmp];
    end
    [Econcat] = solve_l1l2(F,lambda/mu);
    clear F
    start = 1;
    for k=1:K
        E{k} = Econcat(start:start + size(X{k},1) - 1,:);
        start = start + size(X{k},1);
    end
    clear Econcat
    %-------------------4 update G---------------------------------
    Z_tensor = cat(3, Z{:,:});
    Q1_tensor = cat(3, Q1{:,:});
    z = Z_tensor(:);
    clear Z_tensor
    q = Q1_tensor(:);
    [g, objV] = wshrinkObj(z + 1/mu*q,1/mu,sX,0,3);
    G_tensor = reshape(g, sX);
    %-------------------5 update J -------------------------------
    for k=1:K
        D{k}=Z{k}+Q3{k}/mu;
        J{k}=0.5*(D{k}+D{k}');
    end
    %-------------------6 update auxiliary variable---------------
    q = q + mu*(z - g);
    clear g
    clear z
    Q1_tensor = reshape(q, sX);
    clear q
    for k=1:K
        Q2{k} = Q2{k} + mu*(X{k}-X{k}*Z{k}-E{k}); 
        G{k} = G_tensor(:,:,k);
        Q1{k} = Q1_tensor(:,:,k);
        Q3{k} = Q3{k} + mu*(Z{k}-J{k});
    end
    clear G_tensor
    clear Q1_tensor
    %% coverge condition
    Isconverg = 1;
    for k=1:K
        if (norm(X{k}-X{k}*Z{k}-E{k},inf)>epson)
            history.norm_Z = norm(X{k}-X{k}*Z{k}-E{k},inf);
            fprintf('    norm_Z %7.10f    ', history.norm_Z);
            Isconverg = 0;
        end     
        if (norm(Z{k}-G{k},inf)>epson)
            history.norm_Z_G = norm(Z{k}-G{k},inf);
            fprintf('norm_Z_G %7.10f    \n', history.norm_Z_G);
            Isconverg = 0;
        end
    end   
    if (iter>20)
        Isconverg  = 1;
    end
    iter = iter + 1;
    mu = min(mu*pho_mu, max_mu); 
end
S = 0;
for k=1:K
    S = S + abs(Z{k});
end

[~,Ind] = sort(abs(S),1,'descend');
k=MNV(S);
S1 = zeros(N,N);
for i = 1:N
    for r=1:k
        S1(Ind(r,i),i) = S(Ind(r,i),i);
    end
end
W=postprocessor(S1);
for i=1:3
    C1 = SpectralClustering(W,cls_num);
    [ACC(1,i),NMI(1,i),~,~,~,~,~] = AllMeasure(C1,Y);
end
RES=['ACC    ' num2str(mean(ACC)) '    NMI    ' num2str(mean(NMI))];
disp(RES)
