% clc
% clear;

%% 读数据
mpc=case39_UC;
D1=load("load3996.mat");
D=cell2mat(struct2cell(D1))./100;
H=makePTDF(mpc);
ngen=length(mpc.gen(:,1));
%% 决策变量
G=sdpvar(ngen,96);
G1=[zeros(29,96);G];
x=binvar(ngen,96);
cU=sdpvar(ngen,95);

Pgenmax=mpc.gen(:,9)./100;
Pgenmin=mpc.gen(:,10)./100;

Pmax=mpc.branch(:, 6)/100;
Pmin=-Pmax;

qx=1:96;
Qmax=H*D(:,qx)+Pmax;
Qmin=H*D(:,qx)+Pmin;

Mo=mpc.gencost(:,2);
Mc=mpc.gencost(:,3);
%% 目标函数
% cost=mpc.gencost(:,6);
Objective=sum(mpc.gencost(:,6)'*G+mpc.gencost(:,7)'*x)+ones(1,ngen)'*cU*ones(ngen,1);
C=[];
%% 约束
for i=1:96
    for j=1:10
        C=[C,G(j,i)<=x(i,j)*Pgenmax(j),G(j,i)>=x(i,j)*Pgenmin(j)];
    end
    C=[C,H*G1(:,i)<=Qmax(:,i),H*G1(:,i)>=Qmin(:,i),sum(G(:,i))==sum(D(:,i))];
end

for s=2:96
    C=[C,G(:,s)<=0.1*Pgenmax+G(:,s-1),G(:,s)>=-0.1*Pgenmax+G(:,s-1)];
    % for qs=1:10
    %     C=[C,cU(qs,s)>=diff(x,1,2)*Mo(qs),cU(qs,s)>=-diff(x,1,2)*Mc(qs)];
    % end
end
C=[C,cU>=diff(x,1,2).*repmat(Mo,1,95),cU>=-diff(x,1,2).*repmat(Mc,1,95)];
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
    solution=25*value(Objective)+mpc.gencost(:,2)'*gen_stat+24*mpc.gencost(:,7)'*gen_stat;
    
    disp("===============各机组出力情况=================")
    G1=value(G);
    disp(value(G));
    disp("===============总运行成本/美元=================")
    disp(solution);
    disp("===============模型求解速度/s=================")
    disp(etime)
else
    disp('求解出错');
end
sm=1:96;
plot(sm,value(G),LineWidth=1.2);
legend("发电机1","发电机2","发电机3","发电机4","发电机5","发电机6","发电机7","发电机8","发电机9","发电机10")