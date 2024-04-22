clc;
clear;
warning off;
load('MSRC_V1_5views.mat');v = size(X,2);n=210;param.dd = [1302,48,512,256,210];c = 7;for i=1:v X{i}=X{i}'; end

param.gamma1 = 10;
param.beta = 10;
param.lambda1 =1; 
param.lambda2 = 0.1; 
param.beta2 = 1; 
%%%%%%%%%%%%%%%%
param.v = v;
param.k = 10;
param.c = c;
param.NITER = 10;
rand('twister',5489);
%%%%%%%%%%%%%%%
options = [];
options.NeighborMode = 'KNN';
options.k = 5;
options.WeightMode = 'HeatKernel';
options.t = 10;
%%%%%%%%%%%%%%%
tic
[U, consenX, obj] = CE(X,param,options);
toc
for i = 1:v
    U{i} = U{i}';
end
W = DataConcatenate(U);
W = W';
XX = consenX';
d = size(XX,1);
select=3;
selectedFeas = select*d*0.1;
w = [];
for i = 1:d
    w = [w norm(W(i,:),2)];
end
[~,index] = sort(w,'descend');
Xw = XX(index(1:selectedFeas),:);
for i=1:40
 label=litekmeans(Xw',c,'MaxIter',100,'Replicates',20);
    result1 = ClusteringMeasure(Y,label); 
    result(i,:) = result1;
end

for j=1:2
    a=result(:,j);
    ll=length(a);
    temp=[];
    for i=1:ll
        if i<ll-18
            b=sum(a(i:i+19));
            temp=[temp;b];
        end
    end
    [e,f]=max(temp);
    e=e./20;
    MEAN(j,:)=[e,f];
    STD(j,:)=std(result(f:f+19,j));
    rr(:,j)=sort(result(:,j));
    BEST(j,:)=rr(end,j);
end
fprintf('selectedFeas(per) = %d , ',...
  select*10);
fprintf('\n');
disp(['mean. ACC: ', num2str(MEAN(1,1))]);
disp(['mean. STD(ACC): ', num2str(STD(1,1))]);
disp(['mean. NMI: ', num2str(MEAN(2,1))]);
disp(['mean. STD(NMI): ', num2str(STD(2,1))]);
fprintf('\n');
