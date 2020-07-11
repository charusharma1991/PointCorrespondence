function [P,A,ptR,ptC,tim] = matchHasseCone(H1,H2,cliques1,cliques2,coord1,coord2,nodes,nv1,nv2)
node=[];
for i=1:length(nodes)
    n = size(nodes{i},1);
    node = [node n];
end
c1=zeros(sum(node),2);
c2=zeros(sum(node),2);
for i=1:length(node)
    s1=sum(node(1:i-1));
    for j=1:node(i)
        if i<=length(cliques1) && ismember(nodes{i}(j,:),cliques1{i},'rows') && sum(H1(s1+j,:))>0
            [~,id] = ismember(nodes{i}(j,:),cliques1{i},'rows');
            c1(s1+j,:)=coord1{i}(id,:);
        end
        if i<=length(cliques2) && ismember(nodes{i}(j,:),cliques2{i},'rows') && sum(H2(s1+j,:))>0
            [~,id] = ismember(nodes{i}(j,:),cliques2{i},'rows');
            c2(s1+j,:)=coord2{i}(id,:);
        end
    end
end
for i=1:2
    % Calculate the reconstruction matrix
    if i==1
        nonZeroR = find(sum(c1,2));
        d = c1(nonZeroR,:);
        aa = H1(nonZeroR,:);
        aa = aa(:,nonZeroR);
    end
    if i==2
        nonZeroR = find(sum(c2,2));
        d = c2(nonZeroR,:);
        aa = H2(nonZeroR,:);
        aa = aa(:,nonZeroR);
    end
    idsR = setdiff(1:sum(node),nonZeroR);
    ptR{i} = unique(sort(idsR));
    [H{i}, MNeighborhood{i}, MNeighborCoeff{i}] = Conical_Weights( d, aa );
    num{i}=size(H{i},1);
    if ~isempty(ptR{i})
        newH = zeros(sum(node),size(H{i},2));
        k2=1;
        for k1=1:sum(node)
            if ismember(k1,ptR{i})
                newH(k1,:)=0;
            else
                newH(k1,:)=H{i}(k2,:);
                k2=k2+1;
            end
        end
        newH = [newH zeros(sum(node), length(ptR{i}))];
        for k1=1:sum(node)
            if ismember(k1,ptR{i})
                newH=[newH(:,1:k1-1) circshift(newH(:,k1:sum(node)),[0,1])];
            end
        end
        H{i}=newH;
    end
    
end
tic;
for i = 1:1
    for j = 2:2
        for k1 = 1:sum(node)
            for k2 = 1:sum(node)
                C{i}{j}(k1,k2)=norm(H{i}(k1,:)-H{j}(k2,:));
            end
        end
        nonZeroR = find(sum(c1,2));
        nonZeroC = find(sum(c2,2));
        find(sum(coord1{1},2));
        idsR = setdiff(1:sum(node),nonZeroR);%zero rows
        idsC = setdiff(1:sum(node),nonZeroC);%zero cols
        new_cormat = C{i}{j}(nonZeroR,:);
        newnew_cormat= new_cormat(:,nonZeroC);
        c = munkres(newnew_cormat);
        ptR = unique(sort(idsR));
        ptC = unique(sort(idsC));
        if ~isempty(ptR)
            new_c = zeros(1,length(c)+length(ptR));
            k2=1;
            for k1=1:length(new_c)
                if ismember(k1,ptR)
                    new_c(k1)=0;
                    C{i}{j}(k1,:)=4;
                else
                    new_c(k1)=c(k2);
                    k2=k2+1;
                end
            end
            c=new_c;
        end
        if ~isempty(ptC)
            pp=1;
            while(pp<=length(ptC))
                c(c>=ptC(pp)) = c(c>=ptC(pp))+1;
                pp=pp+1;
            end
            for k1=1:sum(node)
                if ismember(k1,ptC)
                    C{i}{j}(:,k1)=4;
                end
            end
        end
        P=c(1:nv1);
        A=c;
        clear c
    end
end
tim=toc;
