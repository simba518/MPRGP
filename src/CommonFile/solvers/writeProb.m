function writeProb(A,b,path);
IA = zeros(size(A,1)+1,1);
JA = zeros(nnz(A),1);
AA = zeros(nnz(A),1);

flag = 0;
counter = 1;
for k=1:size(A,1)
    [IC, JC, CAA] = find(A(k,:));
    if isempty(CAA)
        flag = 1;
    else
        if flag
            IA(k-1)=counter; 
            flag=0;
        end
        IA(k)=counter;
        JA(counter:counter + length(JC)-1)=JC;
        AA(counter:counter + length(JC)-1)=CAA;
        counter = counter + length(JC);
    end
end
IA(end)=IA(1)+nnz(A);

if ~(strcmp('char',class(path)))
    error('2nd argument must be a string that holds a filename');
end    
fd=fopen(path,'w');
if fd<0
	error('Error opening input file!');
end    
fwrite(fd,[size(A,1) nnz(A)],'int32');
for k=1:length(IA)
	fwrite(fd,[IA(k)-1],'int32');
end
for k=1:length(JA)
	fwrite(fd,[JA(k)-1],'int32');
    fwrite(fd,[AA(k)],'double');
end
for k=1:length(b)
	fwrite(fd,[b(k)],'double');
end
fclose(fd);