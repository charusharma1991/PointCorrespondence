function G = hasseDiag(cliques,nodes)
node=[];
for i=1:length(nodes)
    n = size(nodes{i},1);
    node = [node n];
end
%disp(node)
G = zeros(sum(node),sum(node));
% C1 = zeros(sum(node),2);
% C2 = zeros(sum(node),2);
for i = 1:length(node)-1
    %disp(i)
    s1=sum(node(1:i-1));
    for j=1:node(i)
        s2=node(i);
        for k=1:node(i+1)
            if i<length(cliques) && ismember(nodes{i}(j,:),cliques{i},'rows') && ismember(nodes{i+1}(k,:),cliques{i+1},'rows')
                if prod(ismember(nodes{i}(j,:),nodes{i+1}(k,:)))
                    G(s1+j,s1+s2+k)=1;
                    G(s1+s2+k,s1+j)=1;
                end
            end
        end
    end
end
        
end
