function RelativeError = cRerror_Occ(gs,N,n,ptR,ptC)
tau =[];
trueRelativePerm = [];
n_norm = [];
count = 0;
for i = 1:N-1
    for j = i+1:N
    tau = [tau; gs{i}{j}];
    %id1 = setdiff(1:n,find(sum(data{i},2)));
    %id2 = setdiff(1:n,find(sum(data{j},2)));
    %ids = union(id1,id2);
    tru=1:n;
    if ~isempty(ptR{i})
        tru(ptR{i}(ptR{i}<n+1))=0;
    end
    if ~isempty(ptC{j})
        tru(ptC{j}(ptC{j}<n+1))=0;
    end
    %disp(length(tru))
    trueRelativePerm = [trueRelativePerm;tru];
    %a1=setdiff(1:n,ptR{i});
    %a2=setdiff(1:n,ptC{i});
    %count=count+length(union(a1,a2));
    count = count + length(find(tru));
    end
end
%disp()
tau(tau>n)=0;
trueRelativePerm(trueRelativePerm>n)=0;
RelativeError = zeros(size(trueRelativePerm));
ind = find(trueRelativePerm-tau);
RelativeError(ind)=1;
%disp(sum(sum(RelativeError)))
% str = strcat('error_books','.mat');
% save(str)
%RelativeError = sum(sum(RelativeError))/(size(RelativeError,1)*n);
RelativeError = sum(sum(RelativeError))/count;
%save('Omni_err.mat')