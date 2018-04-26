clear
clc

na= 2
nb= 2
nc= 2

kkk= 1;
%do symbolic mesh 
%obtain the topology matrix for A_K

%then assembly the stiffness matrix according to this topology matrix 
%(assume the topology of the stiffness matrix will not change) 

%the sequence of node1 is saved in array 'perm_dof1' [128x1]
%the sequence of node2 is saved in array 'perm_dof2' [128x1]
%the sequence of elements is saved in array 'perm_element' [128x1]
%kpp - whether the next element should be added to the same position 'k' or the next one 'k+1'
 
swch = 1;
IA(1) = 1;		
indx = 2;		 

element = zeros(8,1); %eight elements connected in one node

TFEM = [1,2,3,4,5,6,13,14,15,10,11,12,28,29,30,31,32,33,40,41,42,37,38,29
	4,5,6,7,8,9,16,17,18,13,14,15,31,32,33,34,35,36,43,44,45,40,41,42
	10,11,12,13,14,15,22,23,24,19,20,21,37,38,39,40,41,42,49,50,51,46,47,48
	13,14,15,16,17,18,25,26,27,22,23,24,40,41,42,43,44,45,52,53,54,49,50,51
	28,29,30,31,32,33,40,41,42,37,38,39,55,56,57,58,59,60,67,68,69,64,65,66
	31,32,33,34,35,36,43,44,45,40,41,42,58,59,60,61,62,63,70,71,72,67,68,69
	37,38,39,40,41,42,49,50,51,46,47,48,64,65,66,67,68,69,76,77,78,73,74,75
	40,41,42,43,44,45,52,53,54,49,50,51,67,68,69,70,71,72,79,80,81,76,77,78];

x_e = [ -1 -1 -1 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 -1 -1 -1 
	0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 
	-1 -1 -1 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 -1 -1 -1 
	0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0
	-1 -1 -1 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 -1 -1 -1 
	0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 
	-1 -1 -1 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 -1 -1 -1 
	0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0];

y_e = [ -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 
	-1 -1 -1 -1 -1 -1 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 
	0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 
	0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 
	-1 -1 -1 -1 -1 -1 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 
	-1 -1 -1 -1 -1 -1 0 0 0 0 0 0 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 
	0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 
	0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1];

z_e = [ -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 
	-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0
	-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0
	-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0
	0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1
	0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1
	0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1
	0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1];

perm_element = [];
perm_dof2 = [];
indx_i = [];
indx_j = [];
indx_k = [];
for i=1:81 %nodes
[a_element,b_node] = find(TFEM==i);
A = sortrows([a_element,b_node]);
perm_element = [perm_element;A(:,1)];
perm_dof2 = [perm_dof2;A(:,2)];
indx_i = [indx_i;diag(x_e(A(:,1),A(:,2)))];
indx_j = [indx_i;diag(y_e(A(:,1),A(:,2)))];
indx_k = [indx_i;diag(z_e(A(:,1),A(:,2)))];
end

%element --> node
%1 --> 19
%2 --> 22
%3 --> 16
%4 --> 13
%5 --> 7
%6 --> 10
%7 --> 4
%8 --> 1
				
perm_dof1(find(perm_element==1)) = 19;
perm_dof1(find(perm_element==2)) = 22;
perm_dof1(find(perm_element==3)) = 16;
perm_dof1(find(perm_element==4)) = 13;
perm_dof1(find(perm_element==5)) = 7;
perm_dof1(find(perm_element==6)) = 10;
perm_dof1(find(perm_element==7)) = 4;
perm_dof1(find(perm_element==8)) = 1;
perm_dof1 = perm_dof1';

%find what kkk-th iteration to skip adding kkk, save into RR
k=1
for i=1:191
if perm_element(i)~=perm_element(i+1) & perm_element(i)<perm_element(i+1)
RR(k) = i;
k=k+1;
end
end

k=1
for i=1:191
if perm_element(i)==perm_element(i+1)
RRANY1(k) = i;
k=k+1;
end
end


for k=1:nc
for j=1:nb
for i=1:na

for nxdof=1:3 %next degree of freedom in the same node

for dof=1:192 %total number of dofs for 8 uncoupled elements 

element(:) = -1; %initialize with a non-element value
%-1 indicates that the element does not exist

if (i>1 & j>1 & k>1)
element(1) = (k-2)*(na-1)*(nb-1) + (j-2)*(na-1) + i-1;
end
if (i<na & j>1 & k>1)
element(2) = (k-2)*(na-1)*(nb-1) + (j-2)*(na-1) + i;
end
if (i>1 & j<nb & k>1)
element(3) = (k-2)*(na-1)*(nb-1) + (j-1)*(na-1) + i-1;
end
if (i<na & j<nb & k>1)
element(4) = (k-2)*(na-1)*(nb-1) + (j-1)*(na-1) + i;
end
if (i>1 & j>1 & k<nc)
element(5) = (k-1)*(na-1)*(nb-1) + (j-2)*(na-1) + i-1;
end
if (i<na & j>1 & k<nc)
element(6) = (k-1)*(na-1)*(nb-1) + (j-2)*(na-1) + i;
end
if (i>1 & j<nb & k<nc)
element(7) = (k-1)*(na-1)*(nb-1) + (j-1)*(na-1) + i-1;
end
if (i<na & j<nb & k<nc)
element(8) = (k-1)*(na-1)*(nb-1) + (j-1)*(na-1) + i;
end 

for el=1:8 %eight elements per node
%if (element(el) ~= -1)  %check whether the element exists
if (element(perm_element(dof))~=-1)
	TFEM_CSR(element(perm_element(dof)),perm_dof1(dof)+nxdof-1,perm_dof2(dof)) = 1;%pos %topoloy matrix

	%RR contains nodes after which you do not add kkk but write the DOF into the same position
	if (isempty(find(dof==RR))==1)
	%find the DOF which is added 
	JA(kkk) = 2*((k-1+indx_k(dof))*na*nb + (j-1+indx_j(dof))*na + i+indx_i(dof)) -1*swch
		if swch == 0
			swch = 1
		else
			swch = 0
		end %if
	kkk = kkk+1;
	end %if
end %if 

end %el

end %dof
%next indx meaning next degree of freedom

	IA(indx) = kkk
	indx = indx + 1
end %nxdof

end %i
end %j
end %k


