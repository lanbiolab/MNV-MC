function [W] = postprocessor(Zn)
N = size(Zn,1);
[~,Ind] = sort( Zn,1,'descend' );
for i = 1:N
    Zn(:,i) = Zn(:,i) ./ (Zn(Ind(1,i),i)+eps);
end
%
[uu,s,~] = svd(Zn);
s = diag(s);
r = sum(s>1e-6);
uu = uu(:, 1 : r);
s = diag(s(1 : r));
M = uu * s.^(1/2);
mm = normr(M);
rs = mm * mm';
W = rs.^(2);
end

