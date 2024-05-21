clc
clear;

%% 读数据
mpc=loadcase("case39_ED_updata.m");
load("load3996.mat");
load("Ppvmax.mat");
Ppvmax=Ppvmax/100;
N=10;
T=96;
fpv=0;
Pd=Pd./100;
H=makePTDF(mpc);
Pmax=mpc.branch(:,6)/100;
Pmin=-Pmax;
Qmin=H*Pd+repmat(Pmin,1,T);
Qmax=H*Pd+repmat(Pmax,1,T);
Pgenmax=mpc.gen(:,9)/100;
Pgenmin=mpc.gen(:,10)/100;
nbus=length(mpc.bus(:,1));
ngen=length(mpc.gen(:,1));
cost=mpc.gencost(:,6);
npv=1;
BG=[zeros(nbus-ngen,ngen);eye(ngen)];
Bpv=zeros(nbus,npv);
Bpv(30)=1;

C = [];
Objective = 0;
G=sdpvar(ngen,T,N);
Gpv=sdpvar(npv,T,N);
U=zeros(T,T);
I=eye(T);
v=diag(I);
U=-I+diag(v(1:T-1),-1);

%% 约束
for i=1:N
    % 发电机约束
    C =[C, repmat(Pgenmin,1,T) <= G(:,:,i) <= repmat(Pgenmax,1,T),0 <= Gpv(:,:,i) <= Ppvmax(:,:,i)];
    % 线路潮流约束
    C=[C, Qmin<= H*(BG*G(:,:,i)+Bpv*Gpv(:,:,i)) <= Qmax];
    % 功率平衡约束
    C = [C, sum(G(:,:,i),1)+sum(Gpv(:,:,i),1) == sum(Pd,1)];
    % 爬坡约束
    C = [C,-repmat(0.1*Pgenmax,1,T-1) <= G(:,:,i)*U(:,1:T-1) <= repmat(0.1*Pgenmax,1,T-1)];
    
    Objective = Objective + (1/N)*(sum(cost'*G(:,:,i)+fpv'*Gpv(:,:,i)));
end

% 求解
tic;
reuslt = optimize(C,Objective);
etime=toc;

if reuslt.problem == 0
    
    disp("求解成功");
    gpv=reshape(value(Gpv),T,ngen);
    solution=25*value(Objective);
    disp("===============总运行成本/美元=================")
    disp(solution);
    disp("===============模型求解速度/s=================")
    disp(etime)
    sm=1:96;
    plot(sm,gpv);
else
    disp('求解出错');
end