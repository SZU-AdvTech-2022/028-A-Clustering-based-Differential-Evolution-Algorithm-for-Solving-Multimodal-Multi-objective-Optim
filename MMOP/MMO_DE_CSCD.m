% func_name fname:'MMF1' and so on
% VRmin xl:[1,-1]
% VRmax xu:[3,1]
% n_obj n_obj:2
% NP popsize:200
% Max_Gen Max_Gen:100
function [ps,pf] = MMO_DE_CSCD( func_name,VRmin,VRmax,n_obj,NP,Max_Gen )
ps = cell(Max_Gen,1);
pf = cell(Max_Gen,1);
n_var=size(VRmin,2);    %ȡ���߿ռ��ά�ȣ�ȡVRmin������
Max_FES=Max_Gen*NP;     %�����Ӧ������
%% ��ʼ����Ⱥ
VRmin=repmat(VRmin,NP,1); %% repmat ���ٲ���һ���ϴ�ľ���VRminԭ���Ǿ��߱���������[1,-1]��1��2�У��������200��2��
VRmax=repmat(VRmax,NP,1);
% rand ���������������״Ϊ (200, 2) �൱��200�����壬ÿ�������������������
pos=VRmin+(VRmax-VRmin).*rand(NP,n_var); %initialize the positions of the individuals pos�ǳ�ʼ��Ⱥ
%% ������Ⱥ������ÿ���������Ӧ��fֵ
fitness=zeros(NP,n_obj);   % ��Ӧ��������Ŀ�꺯�������й�
for i=1:NP   %1:200
    fitness(i,:)=feval(func_name,pos(i,:));  % feval����ִ��ָ���ĺ����������ǰ�ֵ���뵽MMF1�����У��õ�f1 f2����Ӧ�ȣ���һ��һ�У�һ������һ�����壩�ر仯
end
fitcount=NP; % ������Ӧ��������Ϊ��Ⱥ��С����200��һ����Ӧ��������Ҫ200��������Ϊһ������������������Ϊ20000/200=100����
G=0; %��ʼ����������

%% ��ѭ��
while fitcount <= Max_FES   % �����ǰ��Ӧ������С�������Ӧ��������200<20000��400<20000
    G=G+1;
    Parentstar=[pos,fitness];  % pos��200��2����һ��������ÿ�����������x�����߱���������Ӧ������y��Ŀ��ֵ��ƴ�ӣ�200��4
    K=10;
    Parentstar = non_domination_scd_kmeans_sort(Parentstar,n_obj,n_var,K); % Parentstar�ǰ�����ӵ���Ⱦ���ֵ�������е����ݣ�n_var+n_obj+1: the level of frontǰ�صȼ�   +2��the decision space distance���߿ռ���� +3��CSCD value���ھ��������ӵ������  +4��the objective space distanceĿ��ռ����
    
    %% ����ÿ��ǰ�صĸ�����
    Total_pf = max(Parentstar(:,n_var+n_obj+1));  % ǰ�صĸ�����������21���ȼ�1-21��
    pop = cell(Total_pf,1);  %pop��21*1��ÿ�ն���һ���վ���ÿ�������ǲ�ͬ�ȼ������и���
    num = zeros(Total_pf,1);  %num��21*1
    
    for i=1:NP % i=1:200
        pop{Parentstar(i,n_obj+n_var+1)} = [ pop{Parentstar(i,n_obj+n_var+1)};Parentstar(i,:)]; %���յȼ�����������ݵ�pop����ȼ�1�����и�������pop��1��1�еľ�������ȼ�2�����и�������pop��2��1�п���
        num(Parentstar(i,n_obj+n_var+1)) = num(Parentstar(i,n_obj+n_var+1)) + 1; %��ʼnum(Parentstar(i,n_obj+n_var+1))��num(1)=0��num(21*1)���ÿ���ȼ��ĸ�����
    end
    
    % ���� CSCD ֵ��ÿ��ǰ���еĸ����������
    for i=1:Total_pf % i=1:21
        [~,index_sorted_on_scd]= sort(pop{i,1}(:,n_var+n_obj+3),'descend');  %pop{1,1}��pop��1��1�У��ȼ�1�������и��壬���յ�7�У�cscdֵ ����ӵ���Ⱦ�����㣩�������򣬷��ض�Ӧ����
        pop{i,1}= pop{i,1}(index_sorted_on_scd,:);  % ����������pop�е�����������򣬷��ذ���7�У�cscdֵ ����ӵ���Ⱦ�����㣩��������Ľ��
    end
    %% DBESM 
    offspring = [];
    for i=1: Total_pf
        u = zeros(num(i),n_var); %num(1)�ǵȼ�1�ĸ�����
        for j= 1 : size(pop{i,1},1) %size(pop{1,1},1)�ǵȼ�1�ĸ�����
            % �����������DBESM����ѡ������
            if i==1  %���ȼ�Ϊ1ʱ
                r_pf(1:3) = 1; %r_pf��1��3�ж�Ϊ1��3*1������[1 1 1]
            else     %���ȼ�����1ʱ
                r_pf = ceil(rand(1,3)*(i-1)); %ceil������ȡ����rand(1,3)����0-1֮�����������i=2ʱ����ȫ����(2-1)����[0 0 0]
            end
            
            if  num(r_pf(1)) < 3 % ����ȼ�1�ĸ�����С��3��i=1��r_pf=[1 1 1]��r_pf(1)��[1 1 1]�еĵ�һ������num(1)�ǵȼ�1�ĸ���������21
                r1 = ceil(rand* num(r_pf(1)));  %rand���Եȼ�1�ĸ�������������ȡ��
            else  % ����ȼ�1�ĸ��������ڵ���3
                r1 = min(ceil(0.1 * num(r_pf(1))),3);  % p=0.1;�ѵȼ�1�ĸ�������С10����Ȼ������ȡ����ȡ����3����Сֵ��һ��r1=3
            end

            % ������룬������Χ���������ھ��߿ռ��ϵľ���
            % [1.79 1.79 1.79;0.67 0.67 0.67](3*2)-[1.79 1.74 1.22;0.67 0.93 -0.79]=[0 0.0508 0.5619;0 -0.2539 1.4701]�����и���ƽ����������ӣ�3*1��
            distance = sum((repmat(pop{i,1}(j,1:n_var), r1 ,1) - pop{r_pf(1),1}(1:r1,1:n_var)).^2,2); %���ȼ�1����ĵ�1������ǰ���У����߱��������Ƴ�r1��3����2�У�pop1��1�У��ȼ�1�����и��壩��ǰ3��ǰ2�У�ǰ��������ľ��߱������������3*2 ��һ��ȫΪ0����ÿ�и���ƽ����������ӣ�3*1��
            if length(distance) == 1 % length(distance)��3
                select_index = r1;
            else
                if length(find(max(distance)==distance))== r1 %find��distance�����ֵ��distance���������length���ж��ٸ���ȵ����ֵ���������3
                    select_index = ceil(rand*r1); %select_index�����������3����rand����һ������������ȡ��
                else
                    distance = max(distance) - distance; %0.6865-[0 0.6865 0.1510](3*1)=[0.6865 0 0.5355]
                    select_pro =  distance./(sum(distance)); %[0.6865 0 0.5355]./(0.6865+0+0.5355)=[0.5618 0 0.4382]
                    true_pro = cumsum(select_pro); %�ۻ��ͣ�[0.5618 0.5618 1]���ڶ�����1 2��ӣ�����������1 2 3���
                    select_index = find(rand<true_pro); %�������С��true_pro�ĸ���������������2 3
                end
            end
            
            r2 = ceil(rand*num(r_pf(2)));  %r_pf(2)=1��num(1)�ǵȼ�1�ĸ����������������������ȡ��
            r3 = ceil(rand*num(r_pf(3)));

            %% ����
            % v��1*2���ȼ�1�ĵ�1������ǰ2��+0.8*���ȼ�1��select_index(1)=2�и���ǰ2��-�ȼ�1��1�и���ǰ2��+�ȼ�1��r2(22)�и���ǰ2��-�ȼ�1��r3(2)�и���ǰ2�У�
            v = pop{i,1}(j,1:n_var) + 0.8*(pop{r_pf(1),1}(select_index(1),1:n_var) - pop{i,1}(j,1:n_var) + pop{r_pf(2),1}(r2,1:n_var) -pop{r_pf(3),1}(r3,1:n_var));
            v = boundConstraint(v, pop{i,1}(j,:), [VRmin(1,:);VRmax(1,:)]); %Խ�������ɲ�Խ�����
            %% ����
            jrand = ceil(rand*n_var); % 2��1
            for h=1:n_var %h=1:2
                if jrand<1 || h==jrand
                    u(j,h) = v(1,h); % �ѱ�����v[1.036 -0.98]�ĵ�һ������ڶ�������u
                else
                    u(j,h) = pop{i,1}(j,h); % u��28*2���ŵȼ�1����1�ĵ�1��x1���2��x2
                end
            end
        end
        offspring = [offspring;u]; % �����Ӵ���u�����ݼ���offspring(2��)�����棬����Ϊǰi���ȼ�������֮�ͣ����϶���200��
    end
    
    % �����Ӵ�
    for i=1:size(offspring,1) % i=1:200
        offspring_fitness(i,:)=feval(func_name,offspring(i,:)); %һ���������߱���x1x2��MMF1������������Ӧ��f�������offspring_fitness����
        fitcount = fitcount +1 ;  % ԭʼfitcount=200��һ�μ������۴�����һ
    end
    temp_pop = [pos,fitness;offspring,offspring_fitness]; %400*4����ų�ʼ��Ⱥ(200*2)����Ӧ��(200*2)���Ӵ���Ⱥ(����һ��)����Ӧ��
    temp_pop = non_domination_scd_kmeans_sort(temp_pop ,n_obj,n_var,K); %400*8���Ӵ�Ҳ������֧��kmeans����
    
    % Environmental seceltion
    Parent = replace_decision_chromosome_kmeans(temp_pop,n_obj,n_var, NP);
    pos = Parent(:,1:n_var); % pos��200*2
    fitness = Parent(:,n_var+1:n_var+n_obj); % Parent���3-4�з�f1��f2
    ps_gen = Parent(Parent(:,n_var+n_obj+1)==1,1:n_var); 
    ps{G,1} = Parent(Parent(:,n_var+n_obj+1)==1,1:n_var); % ȡ��Parent�����5�У��ȼ�������1���У���ȡǰ2�У��ȼ�Ϊ1�ĲŽ�ǰ�ذ���
    pf{G,1} = Parent(Parent(:,n_var+n_obj+1)==1,n_var+1:n_var+n_obj);
end


