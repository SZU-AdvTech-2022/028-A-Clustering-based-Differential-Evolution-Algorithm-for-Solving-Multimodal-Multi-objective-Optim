% Parent = replace_decision_chromosome_kmeans(temp_pop,n_obj,n_var, NP);
% temp_pop intermediate_chromosome �� 400*8���Ӵ�Ҳ������֧��kmeans�����õ���
% NP pop ��200
% M=2��V=2��pop=200
function f = replace_decision_chromosome_kmeans(intermediate_chromosome, M, V, pop)

% ��x1�������У�ѡǰ����
[ C,IA,IC ] = unique(intermediate_chromosome(:,1:V),'rows'); % unique�ǻ�þ����з��ظ����У�����ǰ���е����ݲ���intermediate_chromosome���е�Ψһ�У���û�ظ���������[1 2]��[1 2]���ظ���ֻ��һ�Σ���������1����������C�ǵõ��ĵ�1�����������ݣ�IA�Ƕ�ӦC��һ�е�������IC�Ƕ�Ӧintermediate_chromosome���һ�е�����
intermediate_chromosome = intermediate_chromosome(IA,:);  % �ǰ���1���������������ݣ�400*8��
% IAΪ����C��ȥ���ظ�Ԫ�صľ����е�Ԫ���ھ���A�е�λ�ã���intermediate_chromosome��A��
% ICΪ����A�е�Ԫ���ھ���C�е�λ�ã���C���ļ��У�������������A��
[N, m] = size(intermediate_chromosome); % N=400��m=8
% Get the index for the population sort based on the rank ���ȼ��������У�����
[temp,index] = sort(intermediate_chromosome(:,M + V + 1)); 
clear temp m

% ���ڸ��ݵȼ������Ը���������򣬵õ����ȼ������400�����壨����������������Ӵ���8�����ݣ�
for i = 1 : N
    sorted_chromosome(i,:) = intermediate_chromosome(index(i),:);
end

% Find the maximum rank in the current population �ҳ���ǰ��Ⱥ�е����ȼ�����22
max_rank = max(intermediate_chromosome(:,M + V + 1));

% ��ʼ��ӻ��ڵȼ���ӵ���Ⱦ����ÿ��ǰ�أ�ֱ��������Ⱥ��������
previous_index = 0;
for i = 1 : max_rank
    % Get the index for current rank i.e the last the last element in the
    % sorted_chromosome with rank i. �õ���ǰ�ȼ�Ϊi�ĸ�����
    current_index = max(find(sorted_chromosome(:,M + V + 1) == i)); % find�ҵȼ�Ϊ1�ĸ����������find�õ�81*1����1-81����������sorted_chromosome��ǰ81���ǵȼ�Ϊ1�ĸ��壩��i=2ʱ��current_indexֵΪǰ�����ȼ����������
    
    % ������еȼ�i�ĸ��嶼����ӵ���Ⱥ�У������Ⱥ�Ƿ���䡣
    if current_index > pop  % �����ǰ�ȼ�����������200
        % ����ǣ���ô�ҳ�ǰ�����ȼ����ܸ�����
        remaining = pop - previous_index; % 200��167��ǰ�����ȼ��ĸ������������������������200��������ٸ�����
        temp_pop = ...
            sorted_chromosome(previous_index + 1 : current_index, :);  % �����168-207�У����ĵȼ��ĸ��壩��ֵ��temp_pop��
        [temp_sort,temp_sort_index] = ...
            sort(temp_pop(:, M + V + 3),'descend');  % �ѵ�4�ȼ����壨��ΪҪȡ��4�ȼ���������ĸ��嵽���������200�����壩����7��(����ӵ���Ⱦ���)�������еõ�temp_sort
        for j = 1 : remaining  % j=1:33
            f(previous_index + j,:) = temp_pop(temp_sort_index(j),:); %��f��167�У�ǰ167�а���ǰ�����ȼ��ĸ��壩������ӵ��ĵȼ�������ӵ���Ⱦ��뽵�����е�ǰ33������
        end
        return;
    elseif current_index < pop  % �����ǰ�ȼ�������С��200���������
        % �����еȼ�Ϊi�ĸ���ӵ���Ⱥ�У�һ���ȼ�һ���ȼ��ؼӣ�
        f(previous_index + 1 : current_index, :) = ...
            sorted_chromosome(previous_index + 1 : current_index, :); %һ�ΰѵȼ�i�ĸ���ŵ�f��ۼӣ����ѵ�һ�ȼ�����(1-81��)Ū��f���i=2ʱ���ѵڶ��ȼ�����(82-127��)Ҳ���뵽f�����82-127��
    else
        % �����еȼ�Ϊi�ĸ��������Ⱥ
        f(previous_index + 1 : current_index, :) = ...
            sorted_chromosome(previous_index + 1 : current_index, :);
        return;
    end
    % ��ȡ�����ӵĸ������������1-81ȡ81��82-127ȡ127��
    previous_index = current_index;  % �ѵ�һ�ȼ���������ֵ��previous_index��i=2ʱ��current_indexֵΪǰ�����ȼ����������
end
end


