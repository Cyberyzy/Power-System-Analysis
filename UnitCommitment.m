clc
clear;

%% 读数据
mpc=case39_UC;
load("load3996.mat");
Pd=Pd/100;
H=makePTDF(mpc);
ngen=length(mpc.gen(:,1));
Pgenmax=mpc.gen(:,9)./100;
Pgenmin=mpc.gen(:,10)./100;
Pmax=mpc.branch(:, 6)/100;
Pmin=-Pmax;
T=96;% 时间颗粒度
T0=60;% 启停持续时间
U=zeros(T,T);
U0=zeros(T,T);
I=eye(T);
v=diag(I);
U=-I+diag(v(1:T-1),-1);
for i=1:T0
    U0=U0+diag(v(1:T-i),-i);
end
%% 决策变量
G=sdpvar(ngen,96);
G1=[zeros(29,96);G];
x=binvar(ngen,96);
cU=sdpvar(ngen,95);
Qmin=H*Pd+repmat(Pmin,1,T);
Qmax=H*Pd+repmat(Pmax,1,T);

Mo=mpc.gencost(:,2);
Mc=mpc.gencost(:,3);
%% 目标函数
Objective=sum(100*0.25*mpc.gencost(:,6)'*G)+sum(0.25*mpc.gencost(:,7)'*x)+sum(ones(1,ngen)*cU);
C=[];

%% 约束
%爬坡约束
for s=2:96
    C=[C,G(:,s)<=1.05*Pgenmin+(0.1*Pgenmax-1.05*Pgenmin).*x(:,s-1)+G(:,s-1),G(:,s)>=-1.05*Pgenmin-(0.1*Pgenmax-1.05*Pgenmin).*x(:,s)+G(:,s-1)];
end
% 发电机出力约束
C=[C, repmat(Pgenmin,1,T).*x <= G <= repmat(Pgenmax,1,T).*x];
% 线路潮流约束
C=[C,H*G1<=Qmax,H*G1>=Qmin];
% 有功平衡约束
C=[C, sum(G,1)==sum(Pd,1)];

QS=x*U0;
QL=x*U;
QM=T0*QL+T0;
% 最短启停时间约束
C=[C,cU>=diff(x,1,2).*repmat(Mo,1,95),cU>=-diff(x,1,2).*repmat(Mc,1,95),QM(:,1:T-T0)>=QS(:,1:T-T0)>=T0*QL(:,1:T-T0)];
% 求解
tic;
reuslt = optimize(C,Objective);
etime=toc;
gen_stat=zeros(ngen,1);
for qs=1:ngen
    if(Pgenmax(qs)==0)
        gen_stat(qs)=0;
    else
        gen_stat(qs)=1;
    end
end
if reuslt.problem == 0
    disp("求解成功");
    solution=value(Objective);
    disp("===============总运行成本/美元=================")
    disp(solution);
    disp("===============启停成本/美元=================")
    disp(sum(ones(1,ngen)*value(cU)));
    disp("===============模型求解速度/s=================")
    disp(etime)
else
    disp('求解出错');
end
x1=value(x);
figure
sm=1:96;
plot(sm,value(G),LineWidth=1.2);
xlabel("时间")
ylabel("机组出力情况")
legend("发电机1","发电机2","发电机3","发电机4","发电机5","发电机6","发电机7","发电机8","发电机9","发电机10")
figure
plot(sm,value(x),'LineWidth',1.2);
xlabel("时间")
ylabel("机组启停情况")
legend("发电机1","发电机2","发电机3","发电机4","发电机5","发电机6","发电机7","发电机8","发电机9","发电机10")