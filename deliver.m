clc;
clear;
close all;
%% 题目描述：外卖配送问题
%每一个客户都对应特定的商家
%只有一个骑手 %后续可以换成n个骑手
%无限容量、无限里程、无限时间 %后续可以换成有限
%一般是最小里程等求和 作为目标函数，比较符合自然规律；暂时没想到里程接近
%收到订单后要先到商家才能再到客户

%% 主函数部分

% 数据
Buyer = [-4,-1; -4.5,-0.5; -3,-3; -3,2.5; -3,0.5; -1,-2; -1,-1; -0.5,-1; 
          -0.5,0; 0,3; 1.5,-4; 1,-3; 1, -2; 0.5,-1; 0.5,3; 2.5,-1; 1.5,0; 
          2,1; 3,-1; 3,4];
Seller = [2,-3; 0,4; 1,0; -4,-3; 2, -2; 4,-1; 3,3; 4,2; 4,-1; 3,1;
          -3,-1; 5,0; -2,0; -4,2; -3,1; -3,-2; -1,2; -1,0;0,1;0,2 ];
n = size(Seller,1);


global pre %前提列表，必须要先经过前提
pre(1:n) = 0; %卖家无前提
pre(n+1:2*n) = 1:n; %买家前提是对应卖家
% 测试外卖问题
X = [Seller;Buyer];
[Result,~] = ACO(X, 0);
% 测试TSP
% X = Seller;
% [Result,~] = ACO(X, 1);

%% 画图
plot(Result(:,1),Result(:,2),'o-');
hold on;
plot(Seller(:,1),Seller(:,2),'o');
plot(Result(1,1),Result(1,2),'rp','MarkerSize',9);
%% 尝试添加约束，通解各类TSP、VRP
function [result] = is_available(k,i)
    global Ant_tabu
    global pre
    if ismember(k,Ant_tabu(:,i)) %如果已经路过，则不可行
        result = 0 ;
        return;
    end
    if ismember(pre(k),Ant_tabu(:,i))==0 %如果还没路过商家，则不可行
        result = 0;
        return;
    end
    result = 1;
end

%% 蚁群算法
function [pts, min_distance]=ACO(points,re)
%points是所有点，re表示是否回到起点,0表示不返回
global Ant_tabu
%获取点数量
num = size(points, 1);     %City_Coord的列数即为城市数量
%计算点之间的距离distance矩阵
for i = 1:num
    for j = 1:num
        if i <= j
            distance(i,j) = norm(points(i,:)-points(j,:))+eps;
            distance(j,i) = distance(i,j);
        end
    end
end

%蚂蚁数量
ant_num = 30; 
%初始信息量
t = ones(num, num);
%初始信息量增量
dt = zeros(num, num);
%设置信息量权重,启发量权重,信息素挥发因子,信息素增量常数
a = 1;
b = 5;
p = 0.8;
Q = 200;
%设置最大循环次数
max_iteration = 50;

point_sequence_min=[];

for iteration = 1:max_iteration
    fprintf("aco_iteration : %d\n",iteration);
    %储存蚂蚁走过点的序列
    point_sequence = [];
    %蚂蚁的禁忌列表，表示已经走过的点
    Ant_tabu = zeros(1,ant_num);%0表示广义前驱，放入禁忌列表表示已经经过
    %蚂蚁位置初始化
    point_sequence(1,:) = unidrnd(num/2,1,ant_num);   %随机产生point_sequence的第一行 %肯定是商家
    Ant_tabu(2,:) = point_sequence(1,:);%将初始位置存入禁忌列表中
    %对路径上的信息量作出更新
    t = p * t + dt;                       
    %计算带权重的信息量与启发量的乘积，即概率的分子（numerator）
    numerator = t.^a .* (1./distance).^b;     
    for i = 1:ant_num
        record_distance = 0;
        dt0 = zeros(num, num);
        for j = 1:num-1       
            allow_probability = [];%概率矩阵初始化
            ka = 0;
            for k = 1:num %计算不在禁忌列表中的城市的概率
                ka = ka + 1;
                if  is_available(k,i) %如果k不在蚂蚁i的禁忌列表中
                    allow(ka) = k;
                    allow_probability(ka) = numerator(point_sequence(j,i),k); 
                end
            end
            allow_probability = allow_probability./ sum(sum(allow_probability));  %求出概率
            %轮盘赌法选择下一个点
            allow_probability = cumsum(allow_probability); %列和
            R = rand;
            for kb = 1:ka
                if R < allow_probability(kb)
                    next = allow(kb);
                    break
                end
            end    
            record_distance = record_distance + distance(point_sequence(j,i),next);
            %这边信息素 原来写了距离，笔者认为应该为全局距离的平均值
            %dt0(point_sequence(j,i),next) = Q/distance(point_sequence(j,i),next);%更新从上一个点到下一个点之间路径上的信息量
            point_sequence(j+1,i) = next;
            Ant_tabu(j+2,i) = next;  %将新地点放入禁忌列表
        end
        %该例子中 不返回起点
        if re
            all_record_distance(1,i) = record_distance + distance(point_sequence(num,i),point_sequence(1,i));%计入回到原点的距离
            %dt0(point_sequence(num,i),point_sequence(1,i)) = Q/distance(point_sequence(num,i),point_sequence(1,i));%更新从上一个点回到起点路径上的信息量
            point_sequence(num+1,i) = point_sequence(1,i);%记录到达起点   
        else
            all_record_distance(1,i) = record_distance;
        end
        % 对整条路径添加【平均路径长度】的信息素浓度
        avr_dis = all_record_distance(1,i) / num;
        for j = 1:num-1 
            dt0(point_sequence(j,i),point_sequence(j+1,i)) = Q/avr_dis;
        end
        if re
            dt0(point_sequence(num,i),point_sequence(1,i)) = Q/avr_dis;
        end
        dt = dt + dt0; %累积信息量
         
    end
    [distance_min(iteration), ant_min_i] = min(all_record_distance); %找出走出最短路径的蚂蚁
    point_sequence_min(:,iteration) = point_sequence(:,ant_min_i);%#ok<*AGROW> %存储最短路径
end
[min_distance, min_iteration] = min(distance_min); %#ok<*ASGLU>
final_sequence = point_sequence_min(:,min_iteration);
pts = points(final_sequence,:);
end
