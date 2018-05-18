%clear
%clc

% PARAMETERS TO MODIFY
% --------------------
% DEFINE GRID
%na=2
%nb=2

function [JA,IA,TFEM_CSR,K_diag] = CSR_analysis_2D( na,nb )

% DON'T TOUCH
% -----------
kkk = 1;
swch = 1;
pos = 2;
IA(1) = 1	;
info=0;

% topology matrix of the four elemetns
TFEM = [1,2,3,4,9,10,7,8
3,4,5,6,11,12,9,10
7,8,9,10,15,16,13,14
9,10,11,12,17,18,15,16];

% element coordinates of the four elements
x_e = [-1 -1 0 0 0 0 -1 -1
0 0 1 1 1 1 0 0
-1 -1 0 0 0 0 -1 -1
0 0 1 1 1 1 0 0];

% element coordinates of the four elements
y_e = [-1 -1 -1 -1 0 0 0 0
-1 -1 -1 -1 0 0 0 0
0 0 0 0 1 1 1 1
0 0 0 0 1 1 1 1];

% these vectors show how the property is permuted when going from 1st to end-th node
% how these are related to the node which is being called 
perm_element = [];
perm_dof1 = []; 
perm_dof2 = [];
indx_i = [];
indx_j = [];

for i=1:18 %total number of degrees of freedom for the four elements
[a_element,b_node] = find(TFEM==i);
A = sortrows([a_element,b_node]);
perm_element = [perm_element;A(:,1)];
perm_dof2 = [perm_dof2;A(:,2)];
indx_i = [indx_i;diag(x_e(A(:,1),A(:,2)))];
indx_j = [indx_j;diag(y_e(A(:,1),A(:,2)))];
end

%element --> node (central node)
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

%the number in the brackets are all except these: 	
i=1;
for dof=1:31
	if perm_element(dof)~=perm_element(dof+1)
		if perm_element(dof+1)>=perm_element(dof)
			if indx_i(dof) == indx_i(dof+1) | indx_j(dof) == indx_j(dof+1)
			RR2_excl(i) = dof+1;
			i=i+1
			end
		end
	end
end
RR2 = [1:32];
RR2(RR2_excl)=[];

%-------------------------------------------------------------------------
for j=1:nb
for i=1:na
%element(:) = -1;

for nxdof=1:2
	for dof=1:32
	
		%need it 
		%if any of these numbers is equal to dof then execute the condition: info=0
		%don't count when you switch element 
		%if isempty(find([1,2,3,5,7,8,9,11,13,17,21,23,25,26,27,29,31,32]==dof))==0 %find from only the first DOFs for each kkk
		%info=0;
		%end
		%these numbers are now in RR2

		if isempty(find(RR2==dof))==0
		info=0;
		end

		%if element dof exists then info = 1; 
		%if it does not exist should I add kkk? 
		info_element = 0;
		switch perm_element(dof)
			case 1
				if (i>1 & j>1)
				element = (j-2)*(na-1)+i-1;
				info = 1;
				info_element = 1;
				%else 
				%info = 0;
				end 
			case 2 
				if (i<na & j>1)
				element = (j-2)*(na-1)+i;
				info = 1;
				info_element = 1;
				%else 
				%info = 0;
				end
			case 3 
				if (i>1 & j<nb)
				element = (j-1)*(na-1)+i-1;
				info = 1;
				info_element = 1;
				%else
				%info = 0;
				end
			case 4
				if (i<na & j<nb)
				element = (j-1)*(na-1)+i;
				info = 1;
				info_element = 1;
				%else 
				%info = 0;
				end
				otherwise 
				%info = 0;
		end

% ASSEMBLY THE CSR TOPOLOGY MATRIX
% --------------------------------
if info_element == 1
TFEM_CSR(element, perm_dof1(dof)+nxdof-1, perm_dof2(ndof)) = kkk;
end

% FIND POSITION OF EACH NODE ON THE DIAGONAL OF THE STIFFNESS MATRIX
% ------------------------------------------------------------------
if (nxdof==1 & dof==13)
	K_diag(2*(na*(j-1)+i)-1) = kkk;
end
if (nxdof==2 & dof=17)
	K_diag(2*(na*(j-1)+i)) = kkk;
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
		

	end %dof
	IA(pos) = kkk;
	pos = pos + 1;
end %nxdof

end %na
end %nb

end %function
