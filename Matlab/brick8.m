%brick

B = zeros(6,24);
k_e = zeros(24,24);
C_11 = 10.0;
C_12 = 8.0;
C_44 = 5.0;

a = 2.2/4;
b = 3.12/4;
c = 4.8/4;

%Jacobian = [a 0 0;0 b 0;0 0 c];

%2 integration points
w(1) = 1.0;
w(2) = 1.0;

xi_v(1) = -0.57735026918962576451;
xi_v(2) = -xi_v(1);

eta_v = xi_v;
theta_v = xi_v;

% elastic constants
C = [C_11 C_12 C_12 0 0 0
     C_12 C_11 C_12 0 0 0
     C_12 C_12 C_11 0 0 0
     0 0 0 C_44 0 0 
     0 0 0 0 C_44 0
     0 0 0 0 0 C_44];

for i=1:2 %2 integration points
for j=1:2 %2 integration points
for k=1:2 %2 integration points
    
xi = xi_v(i);
eta = eta_v(j);
theta = theta_v(k);
    
% shape functions:
%N(1) = 1/8*(1-xi).*(1-eta).*(1-theta); 
%N(2) = 1/8*(1+xi).*(1-eta).*(1-theta);
%N(3) = 1/8*(1+xi).*(1+eta).*(1-theta);
%N(4) = 1/8*(1-xi).*(1+eta).*(1-theta);
%N(5) = 1/8*(1-xi).*(1-eta).*(1+theta);
%N(6) = 1/8*(1+xi).*(1-eta).*(1+theta);
%N(7) = 1/8*(1+xi).*(1+eta).*(1+theta);
%N(8) = 1/8*(1-xi).*(1+eta).*(1+theta);

dNdx(1) = -1/8*(1-eta).*(1-theta)/a;
dNdy(1) = -1/8*(1-xi).*(1-theta)/b;
dNdz(1) = -1/8*(1-xi).*(1-eta)/c;

dNdx(2) = 1/8*(1-eta).*(1-theta)/a;
dNdy(2) = -1/8*(1+xi).*(1-theta)/b;
dNdz(2) = -1/8*(1+xi).*(1-eta)/c;

dNdx(3) = 1/8*(1+eta).*(1-theta)/a;
dNdy(3) = 1/8*(1+xi).*(1-theta)/b;
dNdz(3) = -1/8*(1+xi).*(1+eta)/c;

dNdx(4) = -1/8*(1+eta).*(1-theta)/a;
dNdy(4) = 1/8*(1-xi).*(1-theta)/b;
dNdz(4) = -1/8*(1-xi).*(1+eta)/c;

dNdx(5) = -1/8*(1-eta).*(1+theta)/a;
dNdy(5) = -1/8*(1-xi).*(1+theta)/b;
dNdz(5) = 1/8*(1-xi).*(1-eta)/c;

dNdx(6) = 1/8*(1-eta).*(1+theta)/a;
dNdy(6) = -1/8*(1+xi).*(1+theta)/b;
dNdz(6) = 1/8*(1+xi).*(1-eta)/c;

dNdx(7) = 1/8*(1+eta).*(1+theta)/a;
dNdy(7) = 1/8*(1+xi).*(1+theta)/b;
dNdz(7) = 1/8*(1+xi).*(1+eta)/c;

dNdx(8) = -1/8*(1+eta).*(1+theta)/a;
dNdy(8) = 1/8*(1-xi).*(1+theta)/b;
dNdz(8) = 1/8*(1-xi).*(1+eta)/c;

    B(:,:) = 0.0;
    for l=1:8
        B(1,3*l-2) = dNdx(l);
        B(2,3*l-2+1) = dNdy(l);
        B(3,3*l-2+2) = dNdz(l);
        B(4,3*l-2) = dNdy(l);
        B(4,3*l-2+1) = dNdx(l);
        B(5,3*l-2) = dNdz(l);
        B(5,3*l-2+2) = dNdx(l);
        B(6,3*l-2+1) = dNdz(l);
        B(6,3*l-2+2) = dNdy(l);
    end
    
k_e = k_e + B'*C*B*w(i)*w(j)*w(k)*a*b*c;    
    
end
end
end

    