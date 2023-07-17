function [k] = MNV(S)
[~,Ind] = sort(abs(S),1,'descend');
N = size(S,2);
top=floor(N^0.5);
result=zeros(top,1);
SS=zeros(N,top);
for j=2:top
        for i=1:N
            SS(i,j)=find(i==Ind(:,Ind(j,i)));
        end
end

for r=2:top
    tempt=0;  
    for j=2:r
        for i=1:N
            ss=SS(i,j);
            if ss<r
                tempt=tempt+abs(ss-j)/r;  
            else
                tempt=tempt+abs(ss-j)/(N);
            end
        end
        result(r,1)=result(r,1)+tempt*((j/r)^2/r);tempt=0;
    end
    if result(r,1)<result(r-1,1)&&result(r-1,1)<result(r-2,1)
        break;
    end
end
k=r-2;
if r==top
disp('top!')
end
end

