% 5 variantas
%
% 1 lygciu sistema:
% /  2x1 + 5x2 + x3 + 2x4 = 14
% |  -2x1 + 3x3 + 5x4 = 10
% |  x1 - x3 + x4 = 4
% \  -3x1 - 4x2 + x3 + x4 = -6
% Gauso-Zeidelio metodas

% 2 lygciu sistema:
% /  2x1 + 4x2 + 6x3 - 2x4 = 10
% |  x1 + 3x2 + x3 - 3x4 = 0
% |  x1 + x2 + 5x3 + x4 = 10
% \  2x1 + 3x2 - 3x3 - 2x4 = 0
% QR sklaidos metodas





clc, clear all
A=[ 2  5  1  2;
   -2  0  3  5;
    1  0 -1  1;
   -3 -4  1  1]
b=[14;10;4;-6]

n=size(A,1)
Aprad=A;

alpha=[100; 20; 1; 1];  % laisvai parinkti metodo parametrai

Atld=diag(1./diag(A))*A-diag(alpha)
btld=diag(1./diag(A))*b

nitmax=1000;
eps=1e-12;
x=zeros(n,1);x1=zeros(n,1);
fprintf(1,'\n sprendimas iteracijomis:'); 
for it=1:nitmax
    for i=1:n
        x1(i)=(btld(i)-Atld(i,1:i-1)*x1(1:i-1)-Atld(i,i:n)*x(i:n))/alpha(i);
    end
  prec(it)=norm(x1-x)/(norm(x)+norm(x1));
  fprintf(1,'iteracija Nr. %d,  tikslumas  %g\n',it,prec(it))
  if prec(it) < eps, break, end
  x=x1;
end
x
disp('patikrinimas')
Aprad*x-b

% semilogy([1:length(prec)],prec,'r.');grid on,hold on

%--- 1 metodo pabaiga ---

% 2 lygtis (QR-sklaida) 

input('Nor?dami testi, spauskite ENTER')

clc
A=[2  4  6  -2;
   1  3  1  -3;
   1  1  5   1;
   2  3  -3 -2]
Ap=A;  % bus naudojama patikrinimui
b=[10;0;10;0]
n=size(A,1),  nb=size(b,2)

disp(' tiesioginis zingsnis(atspindziai)')
Q=eye(n);
for i=1:n-1
    z=A(i:n,i);
    zp=zeros(size(z));
    zp(1)=sign(z(1))*norm(z);
    omega=(z-zp); omega=omega/norm(omega);
    Qi=eye(n-i+1)-2*omega*omega'
    A(i:n,:)=Qi*A(i:n,:)
    Q=Q*[eye(i-1),zeros(i-1,n-i+1);zeros(n-i+1,i-1),Qi];
end

Q
A

disp('qr skaidos patikrinimas:')
[Qp,Rp]=qr(Ap)


disp('Atvirkstinis zingsnis:')
b1=Q'*b
x=zeros(n,nb);
x(n,:)=b1(n,:)/A(n,n);
for i=n-1:-1:1
    x(i,:)=(b1(i,:)-A(i,i+1:n)*x(i+1:n,:))/A(i,i);
end

x

Ap*x-b
