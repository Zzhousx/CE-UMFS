function [U, ConsenX, obj] = CE(X,param,options)
%%paramter:
% lambda1: local structure preserve
% lambda2: HSIC 
% beta2: nuclear norm
%% ===================== Parameters =====================
gamma1=param.gamma1;
beta=param.beta;
lambda1=param.lambda1;
lambda2=param.lambda2;
beta2=param.beta2;
NITER = param.NITER;
v = param.v;
c = param.c;
%% ===================== initialize =====================
U=cell(1,v);
V=cell(1,v);
L=cell(1,v);
Xnor=cell(1,v);
Xcon=cell(1,v);
E=cell(1,v);
Ei=cell(1,v);
p=cell(1,v);
P=cell(1,v);
qj=cell(1,v);
q=cell(1,v);
Q=cell(1,v);

I=eye(c);
Kn = c;
K_complement = zeros(Kn,Kn);
H = ones(Kn,Kn)*(1/Kn)*(-1) + eye(Kn);
Vcons_length = size(X{1},2);
Vcons = rand(Vcons_length,c);
for vIndex=1:v
    Pn = size(X{vIndex},1);
    P2{vIndex} = eye(Pn);
    P{vIndex} = eye(Pn);
    Q{vIndex} = eye(Pn);
    ni = size(X{vIndex},2);
    U{vIndex}=ones(Pn,c);
    V{vIndex}=ones(ni,c);
    X{vIndex}=NormalizeFea(X{vIndex});
    X{vIndex}=X{vIndex}';
    X{vIndex}=X{vIndex}';
    Xnor{vIndex}=X{vIndex};
    TempNorData=Xnor{vIndex};
    % construct graph 
    Sv = constructW(TempNorData',options);
    D_valuev = sum(Sv,1);
    Dv = spdiags(D_valuev',0,ni,ni);
    Lv = Dv-Sv;
    L{vIndex}=((full(Dv))^(-0.5))*Lv*((full(Dv))^(-0.5)); 
end
for i = 1:v
    Xcon{i} = Xnor{i}';
end
ConsenX = DataConcatenate(Xcon);
%% ===================== updating =====================
for iter = 1:NITER
    for iterv = 1:v
 %% update U
        U{iterv} = U{iterv}.*((P{iterv}*Xnor{iterv}*V{iterv}+P2{iterv}*Xnor{iterv}*Vcons+gamma1*U{iterv})./(P{iterv}*U{iterv}*V{iterv}'*V{iterv}+P2{iterv}*U{iterv}*Vcons'*Vcons+beta*Q{iterv}*U{iterv}+gamma1*U{iterv}*U{iterv}'*U{iterv}+eps));
 %% update V
        % obtain K_complement
        K_complement= K_complement*0;
        for k=1:v
            if (abs(k-iterv)>0) 
            K_complement =  K_complement + H*V{k}'*V{k}*H;
            end
        end
        % update
        V{iterv} = V{iterv}.*((Xnor{iterv}'*P{iterv}*U{iterv})./(V{iterv}*U{iterv}'*P{iterv}*U{iterv}+lambda1*L{iterv}*V{iterv}+lambda2*V{iterv}*K_complement+eps));
 %% construct l_21 norm matrix
        E{iterv} = Xnor{iterv}-U{iterv}*V{iterv}';
        Ei{iterv} = sqrt(sum(E{iterv}.*E{iterv},2)+eps);
        p{iterv} = 0.5./Ei{iterv};
        P{iterv} = diag(p{iterv});
  
        E2{iterv} = Xnor{iterv}-U{iterv}*Vcons';
        Ei2{iterv} = sqrt(sum(E2{iterv}.*E2{iterv},2)+eps);
        p2{iterv} = 0.5./Ei2{iterv};
        P2{iterv} = diag(p2{iterv});
  
        qj{iterv} = sqrt(sum(U{iterv}.*U{iterv},2)+eps);
        q{iterv} = 0.5./qj{iterv};
        Q{iterv} = diag(q{iterv});
    end
 %% update Vcon
    [M,~,N] = svd(Vcons,'econ');
    sum1=zeros(Vcons_length,c);
    sum2=zeros(c,c);
    for sumv=1:v
        sum1=sum1+Xnor{sumv}'*P2{sumv}*U{sumv};
        sum2=sum2+U{sumv}'*P2{sumv}*U{sumv};
    end 
    % update
    Vcons = Vcons.*((sum1)./(Vcons*sum2+beta2*M*N'+eps));
 %% Update gv sum hv sum fv sum

%     gvsum=0;hvsum=0;fvsum=0;jvsum=0;
%     for ghfvIndex=1:v
%         gv{ghfvIndex}=sum(sqrt(sum(U{ghfvIndex}.*U{ghfvIndex},2)));
%         gvsum=gvsum+gv{ghfvIndex}.^(1/(1-r1));
%         
%         hv{ghfvIndex}=trace(V{ghfvIndex}'*L{ghfvIndex}*V{ghfvIndex});
%         hvsum=hvsum+hv{ghfvIndex}.^(1/(1-r2));
%         
%         fv{ghfvIndex}=trace(M'*Vcons*N);
%         fvsum=fvsum+fv{ghfvIndex}.^(1/(1-r3));
%         
%         jv{ghfvIndex}=;
%         jvsum=jvsum+jv{ghfvIndex}.^(1/(1-r4));
%     end
%     
%         fv=trace(M'*Vcons*N);
%         fvsum=fvsum+fv.^(1/(1-r3));
        
%     for abgIndex=1:v
%         gtemp=sum(sqrt(sum(V{abgIndex}.*V{abgIndex},2)));gtemp=gtemp.^(1/(1-r1));
%         beta{abgIndex}=gtemp/gvsum;
        
%         htemp=trace(V{abgIndex}'*L{abgIndex}*V{abgIndex});htemp=htemp.^(1/(1-r2));
%         lambda1{abgIndex}=htemp/hvsum;
        
%         ftemp=trace(U{abgIndex}'*L{abgIndex}*U{abgIndex});ftemp=ftemp.^(1/(1-r3));
%         beta2{abgIndex}=ftemp/fvsum;
%         
%         ftemp=trace(U{abgIndex}'*L{abgIndex}*U{abgIndex});ftemp=ftemp.^(1/(1-r4));
%         beta2{abgIndex}=ftemp/fvsum;
%     end
%         ftemp=trace(M'*Vcons*N);ftemp=ftemp.^(1/(1-r3));
%         beta2=ftemp/fvsum;

%% ===================== calculate obj =====================
     Term5=0;
     for vv=1:v
         K1{vv}=V{vv}'*V{vv};
     end
     for i = 1:v    
        for j = 1:v
                if (abs(j-i)>0)
                    Term5 = Term5 + trace(H*K1{i}*H*K1{j});
                end
        end
     end
     Term6 = beta2*trace(M'*Vcons*N);
    tempobj=0;
    for objIndex=1:v
        Term1 = trace((X{objIndex}-U{objIndex}*V{objIndex}')'*P{objIndex}*(X{objIndex}-U{objIndex}*V{objIndex}'));
        Term2 = beta*sum(sqrt(sum(U{objIndex}.*U{objIndex},2)));
        Term3 = gamma1*trace((U{objIndex}'*U{objIndex}-I)*(U{objIndex}'*U{objIndex}-I)');
        Term4 = lambda1*trace(V{objIndex}'*L{objIndex}*V{objIndex});
        Term7 = trace((X{objIndex}-U{objIndex}*Vcons')'*P2{objIndex}*(X{objIndex}-U{objIndex}*Vcons'));
        tempobj=tempobj+Term1+Term2+Term3+Term4+Term6+Term7;
    end
    tempobj=tempobj+Term6+lambda2*Term5;
    obj(iter)=tempobj;
    if iter == 1
        err = 0;
    else
        err = obj(iter)-obj(iter-1);
    end
    
    fprintf('iteration =  %d:  obj: %.4f; err: %.4f  \n', ...
        iter, obj(iter), err);
    if (abs(err))<1e+0
        if iter > 2
            break;
        end
    end
end










