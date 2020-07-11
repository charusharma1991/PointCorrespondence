function err= GM_RipsSphere_SS(data,N)
k=5;
P={};
for ii=1:1
for m=1:N-1
    for n=m+1:N
        data1=data{m};
        data2=data{n};
        nv1 = length(data1);
        nv2 = length(data2);
        [G1,G2] = RipsSphere(data1,data2,nv1,nv2,k);
        [~,cl1,~] = k_clique(4,G1);
        [~,cl2,~] = k_clique(4,G2);

        cliques1{m}{1}=find(sum(data1,2));
        cliques2{n}{1}=find(sum(data2,2));
        coord1{m}{1}=data1(cliques1{m}{1},:);
        coord2{n}{1}=data2(cliques2{n}{1},:);
        count=1;
        for i=1:nv1
            for j=i+1:nv1
                if G1(i,j)==1
                    cliques1{m}{2}(count,:)=[i,j];
                    coord1{m}{2}(count,:)=sum([data1(i,:);data1(j,:)])/2;
                    count=count+1;
                end
            end
        end
        
        count=1;
        for i=1:nv2
            for j=i+1:nv2
                if G2(i,j)==1
                    cliques2{n}{2}(count,:)=[i,j];
                    coord2{n}{2}(count,:)=sum([data2(i,:);data2(j,:)])/2;
                    count=count+1;
                end
            end
        end
        count1=1;
        for i=1:length(cl1)
            if length(cl1{i})==3
                cliques1{m}{3}(count1,:)=cl1{i};
                coord1{m}{3}(count1,:) = sum(data1(cl1{i},:))/3;
                count1=count1+1;
            end
        end
        count1=1;
        for i=1:length(cl2)
            if length(cl2{i})==3
                cliques2{n}{3}(count1,:)=cl2{i};
                coord2{n}{3}(count1,:) = sum(data2(cl2{i},:))/3;
                count1=count1+1;
            end
        end
        nodes{m}{n}{1}=(1:nv1)';
        for i=2:3
            if (i<=length(cliques1{m}) && i<=length(cliques2{n}))
                nodes{m}{n}{i} = union(cliques1{m}{i},cliques2{n}{i},'rows');
            else
                if length(cliques1{m})>length(cliques2{n}) && i>length(cliques2{n})
                    nodes{m}{n}{i}=cliques1{m}{i};
                elseif length(cliques1{m})<length(cliques2{n}) && i>length(cliques1{m})
                    nodes{m}{n}{i}=cliques2{n}{i};
                end
            end
        end
        H1{m} = hasseDiag(cliques1{m},nodes{m}{n});
        H2{n} = hasseDiag(cliques2{n},nodes{m}{n});
        [P{m}{n},A{m}{n},ptR{m},ptC{n},tim]=matchHasseSphere(H1{m},H2{n},cliques1{m},cliques2{n},coord1{m},coord2{n},nodes{m}{n},nv1,nv2);
        clear cliques1 cliques2 coord1 coord2 data1 data2 nv2 G1 G2 cl1 cl2 nodes H1 H2 count count1
    end
end
err = cRerror_Occ(P,N,nv1,ptR,ptC);
end
