na=2
nb=2

element = zeros(4,1);
kkk = 1;
swch = 1;



TFEM = [1,2,3,4,9,10,7,8
        3,4,5,6,11,12,9,10
        7,8,9,10,15,16,13,14
        9,10,11,12,17,18,15,16];
        
x_e = [-1 -1 0 0 0 0 -1 -1
       0 0 1 1 1 1 0 0 
       -1 -1 0 0 0 0 -1 -1
       0 0 1 1 1 1 0 0];
       
y_e = [-1 -1 -1 -1 0 0 0 0 
       -1 -1 -1 -1 0 0 0 0 
       0 0 0 0 1 1 1 1 
       0 0 0 0 1 1 1 1];
       
perm_element = [];
perm_dof2 = [];
indx_i = [];
indx_j = [];

for i=1:18 %nodes
[a_element,b_node] = find(TFEM==i);
A = sortrows([a_element,b_node]);
perm_element = [perm_element;A(:,1)];
perm_dof2 = [perm_dof2;A(:,2)];
indx_i = [indx_i;diag(x_e(A(:,1),A(:,2)))];
indx_j = [indx_j;diag(y_e(A(:,1),A(:,2)))];
end

%element --> node
%1 --> 5
%2 --> 7
%3 --> 3
%4 --> 1
				
perm_dof1(find(perm_element==1)) = 5;
perm_dof1(find(perm_element==2)) = 7;
perm_dof1(find(perm_element==3)) = 3;
perm_dof1(find(perm_element==4)) = 1;
perm_dof1 = perm_dof1';

%find what kkk-th iteration to skip adding kkk, save into RR
k=1
for i=1:31
if perm_element(i)~=perm_element(i+1) & perm_element(i)<perm_element(i+1)
RR(k) = i;
k=k+1;
end
end

%-------------------------------------------------------------------------

for j=1:nb
for i=1:na

element(:) = -1;
 if (i>1 & j>1)
 	element(1) = (j-2)*(na-1) + i-1;
 end
 if (i<na & j>1)
	element(2) = (j-2)*(na-1) + i;
 end
 if (i>1 & j<nb)
 	element(3) = (j-1)*(na-1) + i-1;
 end
 if (i<na & j<nb)
	element(4) = (j-1)*(na-1) + i;
 end

for nxdof=1:2
for dof=1:32

if isempty(find([1,2,3,5,7,8,9,11,13,17,21,23,25,26,27,29,31,32]==dof))==0
      info=0;
end

%for el=1:4
	if (element(perm_element(dof))~=-1)
		info = 1;	%element does exist
	else
	 	info = 0;	%element does not exist
	end

		if (isempty(find(dof==RR))==1 & info==1)		
			%find the DOF which is added 
				JA(kkk) = 2*((j-1+indx_j(dof))*na + i+indx_i(dof)) -1*swch
			if swch == 0
				swch = 1;
			else
				swch = 0;
			end %if
		kkk = kkk+1;
		end %if
	
%end %el


end %dof
end %nxdof

end %na
end %nb
