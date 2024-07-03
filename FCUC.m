%% 全局信息

clc
clear;
mpc = loadcase('case39_FCUC');
load('load3996.mat');
load('Ppvmax.mat');
locat_pv = [1 30;2 31;3 32];%光伏所在节点
num_pv = size(locat_pv,1);
load = Pd/mpc.baseMVA;%负载标幺化;
Period_num = size(load,2);%时间段数

%% mpc相关数据

Gen_num = size(mpc.gen,1);%发电机台数
Node_num = size(mpc.bus,1);%节点数
Gen_bus = [(1:Gen_num).' mpc.gen(:,1)];%发电机所在节点
Gen_enabled = double(mpc.gen(:,9)~=0);%未停用的发电机
Pmax = mpc.gen(:,9);
Qmax = mpc.gen(:,4);
S = Pmax;%机组最大出力

Gen_site = zeros(Node_num,Gen_num);%Gen_site矩阵可以对g向量实现补零
for j = 1:Gen_num
    Gen_site(Gen_bus(j,2), Gen_bus(j,1)) = 1;
end

Gen_max = mpc.gen(:,9)/mpc.baseMVA;
Gen_max = repmat(Gen_max, 1, Period_num);%一维变二维，为了能直接与g_Gen比较
Gen_min = mpc.gen(:,10)/mpc.baseMVA;
Gen_min = repmat(Gen_min, 1, Period_num);%发电机出力限制
Branch_max = mpc.branch(:,6)/mpc.baseMVA;
Branch_max = repmat(Branch_max, 1, Period_num);
Branch_min = -Branch_max;%支路潮流限制
H = makePTDF(mpc);%功率传输分布因子矩阵
c = mpc.gencost(:,6);%成本系数（一次项）
cost_const = mpc.gencost(:,7);%成本系数（常数项）
start_cost = mpc.gencost(:,2);%启动成本
shut_cost = mpc.gencost(:,3);%关停成本

%% 频率安全相关数据

D  = [2.7 2.9 2.1 2.5 2.7 2.6 2.4 2.8 2.2 2.6].';
KG = [15  16  22  21  20  28  23  17  18  21].';
Tj = [15  14  13  13  15  15  10  9   9   15].';
FR = [0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2].';
RoCoFm = 0.5/50;%最大频率变化率,Hz/s
fss = 0.3/50;%最大稳态频率偏差,Hz
fcr = -0.5/50;%最大频率偏差,Hz
dPL = 250/mpc.baseMVA;

TR = 12;

alpha = (pi * fcr * S .* Tj) ./ (4 * dPL * TR * mpc.baseMVA);

%% 光伏相关数据

Ppv_max = reshape(Ppvmax,96,1000);
Ppv_max = Ppv_max/mpc.baseMVA;%光伏标幺化
Ppv_site = zeros(Node_num,num_pv);%光伏位于3个节点上面
for j = 1:num_pv
    Ppv_site(locat_pv(j,2), locat_pv(j,1)) = 1;
end
[idx,CP]=kmeans(Ppv_max',Gen_num);
CPs=reshape(CP',1,96,Gen_num);
% epsilon = 0.05;
% z = norminv(1-epsilon, 0, 1);
% N = 50;
% [idx,Ppv_max] = kmeans(Ppv_max.', N);
% u = mean(Ppv_max, 1);
% sigma = sqrt(var(Ppv_max, 0, 1));

% Ppv_max = 1 * repmat(u - sigma * z, num_pv, 1);%增大上限，出现弃光
punishment_cost = 3750;%弃光惩罚系数,(0.15*1000*100*0.25),单位：$/(100MW*0.25h)

%% 机组启停相关数据

Ru =  0.1*Gen_max(:,  2:Period_num);
Rd = -0.1*Gen_max(:,  1:Period_num-1);
Su =  1.05*Gen_min(:, 2:Period_num);
Sd = -1.05*Gen_min(:, 1:Period_num-1);
To = 4;
Tc = 4;
Uo = zeros(Period_num,Period_num-To);
for i=1:Period_num
    for j=1:Period_num-To
        if (i > j)&&(i-j <= To)
            Uo(i,j)=1;
        end
    end
end
Uc = zeros(Period_num,Period_num-Tc);
for i=1:Period_num
    for j=1:Period_num-Tc
        if (i > j)&&(i-j <= Tc)
            Uc(i,j)=1;
        end
    end
end

%% 构造约束条件

%决策变量
g_Gen =  sdpvar(Gen_num,  Period_num,  'full');
P_node = sdpvar(Node_num, Period_num,  'full');
c_oc  =  sdpvar(Gen_num,  Period_num-1,'full');
x_Gen =  binvar(Gen_num,  Period_num,  'full');
z_Gen =  diff(x_Gen, 1, 2);
Ppv = sdpvar(num_pv,Period_num,'full');

%约束条件
Constraints = [];
%功率平衡
Constraints = [Constraints, Gen_site * g_Gen + Ppv_site * Ppv - load == P_node];
Constraints = [Constraints, sum(P_node) == 0];
%发电机出力限制
Constraints = [Constraints, x_Gen .* Gen_min <= g_Gen <= x_Gen .* Gen_max];
%支路潮流限制
Constraints = [Constraints, Branch_min <= H * P_node <= Branch_max];
%爬坡限制
Constraints = [Constraints, diff(g_Gen,1,2) <= Su + (Ru-Su) .* x_Gen(:, 1:Period_num-1)];
Constraints = [Constraints, diff(g_Gen,1,2) >= Sd + (Rd-Sd) .* x_Gen(:, 2:Period_num)];
%启停成本约束
Constraints = [Constraints, c_oc >=  z_Gen .* repmat(start_cost, 1, Period_num-1)];
Constraints = [Constraints, c_oc >= -z_Gen .* repmat(shut_cost, 1, Period_num-1)];
%机组最小持续开停机时间约束
Constraints = [Constraints, x_Gen * Uo       >=  To * z_Gen(:, 1:Period_num-To)];
Constraints = [Constraints, (1 - x_Gen) * Uc >= -Tc * z_Gen(:, 1:Period_num-Tc)];
Constraints = [Constraints, z_Gen(: , Period_num-To+1 :Period_num-1) == 0];
Constraints = [Constraints, z_Gen(: , Period_num-Tc+1 :Period_num-1) == 0];
%光伏出力限制
% Constraints = [Constraints, 0 <= Ppv <= Ppv_max];
Constraints = [Constraints, 0 <= Ppv <= 4*CP(1:3,:)];
%频率最低点约束
S1 = repmat(S, 1, Period_num);
Ssys = mpc.baseMVA;
Tj_sys = (x_Gen .* S1)' * Tj ./Ssys;
D_sys  = (x_Gen .* S1)' * D  ./Ssys;
KG_sys = (x_Gen .* S1)' * KG ./Ssys;

dP_Gm = repmat(KG * abs(fcr), 1, Period_num) - ...
    ((KG * abs(fcr)) .* (ones(Gen_num,1) - FR)) * ...
    (ones(1,Period_num) - alpha.' * x_Gen);
dPG_r = min(dP_Gm, Gen_max .* x_Gen - g_Gen);
Constraints = [Constraints, sum(dPG_r, 1) >= repmat(dPL, 1, Period_num) + 1.2*(D_sys.') * fcr];
%频率变化率约束
Constraints = [Constraints, -RoCoFm * Tj_sys.' <= repmat(dPL, 1, Period_num) <= RoCoFm * Tj_sys.'];
% 准稳态频率约束
Constraints = [Constraints,-fss * (KG_sys + D_sys).' <= repmat(dPL, 1, Period_num) <= fss * (KG_sys + D_sys).'];

%目标函数
P_cost = sum(c.'*g_Gen)*mpc.baseMVA*0.25 + sum(cost_const.'*x_Gen)*0.25 + sum(c_oc(:));
%弃光惩罚函数
P_cost = P_cost + 3750 * sum(4*CP(1:3,:) - Ppv, "all");
%% 使用gurobi求解

%设置求解器
ops = sdpsettings('verbose', 1,'debug',1);
% ops = sdpsettings('verbose', 0, 'solver', 'gurobi');
%求解
solve = optimize(Constraints, P_cost, ops);
result = value(g_Gen);
solvertime = solve.solvertime;
total_cost = value(P_cost);
start_shut_cost = sum(value(c_oc), "all");
Ppv = value(Ppv);
curtailment = CP(1:3,:) - Ppv;

%% 结果输出

disp("========优化求解用时(s)：========");
disp(solvertime);
disp("========当日总成本($)：========");
disp(total_cost);
disp("========机组启停成本($)：========");
disp(start_shut_cost);
t = 0.25:0.25:24;
subplot(2,1,1);
plot(t, mpc.baseMVA * result, LineWidth=1.2);
legend("机组1","机组2","机组3","机组4","机组5","机组6","机组7","机组8","机组9","机组10");
title("24小时内各机组出力变化曲线");
xlabel("时间/h");
ylabel("有功出力/MW");
subplot(2,1,2);
plot(t, mpc.baseMVA * [Ppv; curtailment] ,LineWidth=1.2);
legend("光伏30","光伏31","光伏32","弃光30","弃光31","弃光32");
title("24小时内各机组出力变化曲线");
xlabel("时间/h");
ylabel("有功出力/MW");