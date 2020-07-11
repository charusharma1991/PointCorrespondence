function err= GM_RipsEllipse_SS(data,N)
%Input: data cell of size N contains each matrix of (x,y) coordinates, N is number of sets (images) to match
%output: err: error percentage of matching
k=5;
P={};
for ii=1:1
for m=1:N-1
    for n=m+1:N
        data1=data{m};
        data2=data{n};
        nv1 = length(data1);
        nv2 = length(data2);
        [G1,G2] = RipsEllipse(data1,data2,nv1,nv2,k);
        [~,cl1,~] = k_clique(4,G1);
        [~,cl2,~] = k_clique(4,G2);

        cliques1{1}=find(sum(data1,2));
        cliques2{1}=find(sum(data2,2));
        coord1{1}=data1(cliques1{1},:);
        coord2{1}=data2(cliques2{1},:);
        count=1;
        for i=1:nv1
            for j=i+1:nv1
                if G1(i,j)==1
                    cliques1{2}(count,:)=[i,j];
                    coord1{2}(count,:)=sum([data1(i,:);data1(j,:)])/2;
                    count=count+1;
                end
            end
        end
        
        count=1;
        for i=1:nv2
            for j=i+1:nv2
                if G2(i,j)==1
                    cliques2{2}(count,:)=[i,j];
                    coord2{2}(count,:)=sum([data2(i,:);data2(j,:)])/2;
                    count=count+1;
                end
            end
        end
        count1=1;
        for i=1:length(cl1)
            if length(cl1{i})==3
                cliques1{3}(count1,:)=cl1{i};
                coord1{3}(count1,:) = sum(data1(cl1{i},:))/3;
                count1=count1+1;
            end
        end
        count1=1;
        for i=1:length(cl2)
            if length(cl2{i})==3
                cliques2{3}(count1,:)=cl2{i};
                coord2{3}(count1,:) = sum(data2(cl2{i},:))/3;
                count1=count1+1;
            end
        end
        nodes={};
        nodes{1}=(1:nv1)';
        for i=2:3
            if (i<=length(cliques1) && i<=length(cliques2))
                nodes{i} = union(cliques1{i},cliques2{i},'rows');
            else
                if length(cliques1)>length(cliques2) && i>length(cliques2)
                    nodes{i}=cliques1{i};
                elseif length(cliques1)<length(cliques2) && i>length(cliques1)
                    nodes{i}=cliques2{i};
                end
            end
        end
        H1 = hasseDiag(cliques1,nodes);
        H2 = hasseDiag(cliques2,nodes);
        [P{m}{n},A{m}{n},ptR{m},ptC{n},tim]=matchHasseEllipse(H1,H2,cliques1,cliques2,coord1,coord2,nodes,nv1,nv2);
        clear cliques1 cliques2 coord1 coord2 data1 data2 nv2 G1 G2 cl1 cl2 nodes H1 H2 count count1
    end
end
err = cRerror_Occ(P,N,nv1,ptR,ptC);
end
