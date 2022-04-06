%
clear;
board_size=[10,4];
mul=board_size(1)*board_size(2);
num_group=mul/20;
h=board_size(1)/2;
w=board_size(2)/2;
c_shape=cell(mul/4,1);
for i = 1:h
    for j=1:w
        c_shape{(i-1)*w+j}=[1,1;2,1;1,2;2,2]+[2*i-2,2*j-2];
    end
end
        

board_map=shapes2map(c_shape,board_size);
rank0=calFitness(c_shape);
rank=rank0;

num_pop=5;
for i=1:num_pop
    v_c_shape{i}=c_shape;
end


n=0;
while rank~=0
    n=n+1;
    fprintf('n: %d\t',n);
    v_rank=zeros(num_pop,1);
    for i=1:num_pop
       v_rank(i)=calFitness(v_c_shape{i});
    end
    v_p=exp(-v_rank/rank0*2);%test
    v_p_1=v_p/sum(v_p);
    cum_v_p_1=cumsum(v_p_1);
    rn=rand;
    index_pop=0;
    for i=1:num_pop
       if rn<=cum_v_p_1(i)
           index_pop=i;
           break;
       end
    end
    c_shape=v_c_shape{index_pop};
    
%     small_map=shape2SmallMap(c_shape{1});
%     name=getSmallMapName(small_map);
    [s1,s2,shape_index]=selectShape(c_shape,board_map);
    [new_s1,new_s2]=crossoverShape(s1,s2);
    % update c_shape and board_map
    c_shape{shape_index(1)}=new_s1;
    c_shape{shape_index(2)}=new_s2;
    board_map=shapes2map(c_shape,board_size);
    rank=calFitness(c_shape);
    fprintf('rank: %d',rank);
    
    %update population
    rn2=rand;
    [max_rank,max_index]=max(v_rank);
    if rank<=max_rank ||rn2<0.5
        flag=1;
        for i=1:num_pop
            if isequal(v_c_shape{i},c_shape)
                flag=0;
                disp('  equal');
                break;
            end
        end
        if flag==1
            
%             v_c_shape{max_index}=c_shape;
            v_p2=exp(v_rank/rank0/2);%test
            v_p_1_2=v_p2/sum(v_p2);
            cum_v_p_2=cumsum(v_p_1_2);
            rn3=rand;
            index_pop2=0;
            for i=1:num_pop
               if rn3<=cum_v_p_2(i)
                   index_pop2=i;
                   break;
               end
            end
            v_c_shape{index_pop2}=c_shape;
            fprintf('\tupdate pop %d',index_pop2);
        end
    end
    fprintf('\n');
end
figure();
image(board_map,'CDataMapping','scaled')
axis equal;
%% GA functions
% selection
function [s1,s2,shape_index]=selectShape(c_shape,board_map)
    % if the num of the shape more than the target num, the chance is larger
    
    %
%     shape_index=[0,0];
    p2=[0,0];
    while true
        rn=randi(length(c_shape));
        s1=c_shape{rn};

        % should only select the shapes next to each other
        rn2=randi(4);% one of the 4 position of the shape
        p1=s1(rn2,:);
        rn3=randi(4);% one of the 4 direction
        if rn3==1
            p2=p1+[1,0];
        elseif rn3==2
            p2=p1-[1,0];
        elseif rn3==3
            p2=p1+[0,1];
        elseif rn3==4
            p2=p1-[0,1];
        end
        % check if p2 is inside the board
        if min(p2)>=1 && p2(1)<=size(board_map,1) && p2(2)<=size(board_map,2)
            % check if p2 belong to another shape
            if board_map(p1(1),p1(2))~=board_map(p2(1),p2(2))
                break;
            end
        end
    end
    s2=c_shape{board_map(p2(1),p2(2))};
    shape_index=[rn,board_map(p2(1),p2(2))];
    
end


% small crossover
function [new_s1,new_s2]=crossoverShape(s1,s2)
    m_visited=zeros(4,4);
    while true
        rn = randi([1 4],1,2);% make try more than 1 position at a time.
        % check if visied before
        if m_visited(rn(1),rn(2))==1
            continue;
        end
        m_visited(rn(1),rn(2))=1;
        new_s1=s1;
        new_s2=s2;
        % swap p between s1 and s2
        p1=s1(rn(1),:);
        p2=s2(rn(2),:);
        new_s1(rn(1),:)=p2;
        new_s2(rn(2),:)=p1;
        
        % check if new s1 and new s2 are basic shapes
        if ~isempty(getShapeName(new_s1))&& ~isempty(getShapeName(new_s2))
            break;
        end
        % check if all possible cross is tried
        if isequal(m_visited,ones(4,4))
            % don't change
            new_s1=s1;
            new_s2=s2;
            break;
        end
    end
end

% big crossover among several shapes
% function

% evaluate fitness
function rank=calFitness(c_shape)
    % cal the num of the basic shapes
    vname=zeros(1,length(c_shape));
    for i =1:length(c_shape)
        name=getShapeName(c_shape{i});
        vname(i)=name;
    end
    vnum=zeros(1,5);
    vnum(1)=length(strfind(vname,'T'));
    vnum(2)=length(strfind(vname,'I'));
    vnum(3)=length(strfind(vname,'O'));
    vnum(4)=length(strfind(vname,'L'));
    vnum(5)=length(strfind(vname,'Z'));
    target_num=length(c_shape)/5;
    rank=sum(abs(vnum-target_num));
    
%     %for test
%     rank=abs(vnum(4)-2)+abs(vnum(3)-1);
end
%% basic functions

function F=getSmallMapName(small_map)
    F='';
    vletter=['T','I','O','L','Z'];
    find_flag=-1;
    for i =1:length(vletter)
        l=vletter(i);
        c_smallMap=getSmallMap(l);
        for j=1:length(c_smallMap)
            find_flag=isequal(small_map,c_smallMap{j});
            if find_flag==1
                F=l;
                break;
            end
        end
        if find_flag==1
            break;
        end
    end
end

function F=shape2SmallMap(s)
%     s=c_shape{1};
    s=s-min(s)+[1,1];
    xy=max(s);
    small_map=zeros(xy);
    for j =1:size(s,1)
        p=s(j,:);
        small_map(p(1),p(2))=1;
    end
    F=small_map;
end

%define shape
function F=getSmallMap(flag)
    if flag=='T'
        temp=ones(2,3);
        temp(2,1)=0;
        temp(2,3)=0;
        temp2=flip(temp,1); 
        F{1}=temp;
        F{2}=temp2';
        F{3}=temp2;
        F{4}=F{1}';
    elseif flag=='I'
        temp=ones(1,4);
        F{1}=temp;
        F{2}=temp';
    elseif flag=='O'
        F{1}=ones(2,2);
    elseif flag=='L'
        temp=ones(3,2);
        temp(1:2,2)=zeros(2,1);
        temp2=flip(temp,1); %t-d
        temp3=flip(temp,2); %l-r
        temp4=flip(temp3,1);
        F{1}=temp;
        F{2}=temp2';
        F{7}=temp2;
        F{8}=F{1}';
        F{5}=temp3;
        F{6}=temp4';
        F{3}=temp4;
        F{4}=F{5}';
    elseif flag=='Z'
        temp=ones(3,2);
        temp(1,1)=0;
        temp(3,2)=0;
        temp2=flip(temp,1); 
        F{1}=temp;
        F{2}=temp2';
        F{3}=temp2;
        F{4}=F{1}';
    end
end

function F=getShapeName(s)
    F=getSmallMapName(shape2SmallMap(s));
end

function F=shapes2map(c_shape,board_size)
    board_map=zeros(board_size);
    for i=1:length(c_shape)
        s=c_shape{i};
        for j =1:size(s,1)
            p=s(j,:);
            board_map(p(1),p(2))=i;
        end
    end
    F=board_map;
end