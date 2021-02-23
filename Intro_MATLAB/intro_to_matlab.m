%click Run Section to run each demo section
%Or use the keyboard shortcut
%Or type them in command window(recommanded)

%%
%define a variable and basic commands (page 6)
clear all ;clc %this clears all existing variables and cleans the screen
n=4;
m=2;
k=n^m
%%

%%
clc
var1=rand

%%
clc
format long  %display more digits
var1
var1short=single(var1) %convert to single precision
format short
%%
clc
m=m+1  %reassign a variable
%m=M+1 %matlab is case sensitive!!!

astring='today the sky is blue';  %matlab takes integer, real, char, etc. no need to define

pwd %prints your current path
ls
%%
help sum  %short manual
doc sum   %long manual
lookfor variance  %search for

%%
clear all;clc
%something about vectors (page 7)
t1=[1 2 3 4 5]; %row vector
t2=[1;2;3;4;5]; %column vector
%%
t3=1:5; %implicit loop, default det1lta is 1
length(t1)
t4=1:0.5:5;
length(t4)

t5=linspace(1,5,10);

c=rand(1,10);

c_new1=c(1:3);
c_new2=c(3:end); %=c(3:length(c))
c_new3=c(:);

c' %transpose the vector
a=c;
a=a*3
%%
clc
anew=a.*c
try
    B=a*c %this is wrong
catch error
    error
    pause
end

B2=a*c'
B3=a'*c

%%
%input matrix manually  (page 8)
clear all; clc;
a=[1 2 3;4 5 6;7 8 9;10 11 12]

%%
clc
a=zeros(3);
b=ones(3);
c=rand(3);
d=eye(3);
e=diag([1 2 3])

abcd=[a b;c d]

size(abcd)

a(1)=5;  %when using a sinlge index o
a(2)=10
a(1,3)=7;
a(3,3)=6


rr=a(2,:)
a(:,3)
a(:,1:3)
a(1:2,1:2)


size(rr)
[m,n]=size(rr)


%%
% element-wise operation
clear all; clc
A=ones(4)
B=rand(4)
C1=A+A
C2=A*3
C3=A.*A
C4=A*A
C5=A./B
D=sin(A)

%%
%matrix operation
clear all; clc
A=eye(3)
B=rand(3)
C1=A*A
C2=A*B
C3=A/B
C4=A*inv(B)
max(max(C3-C4))
%%
clc
A=eye(1000);
B=rand(1000);
C3=A/B;
C4=A*inv(B);
max(max(C3-C4))
%%
clc
B=rand(5);
max(B)
B   %change direction  row / column
min(B)
min(B,[],2)

sum(B)
sum(B,2)

diag(B)
%%
clc
B
ind=find(B>0.5)
[m,n]=find(B>0.5)
T=B>0.5  %return a matrix of 1 and 0
B2=B.*T %extract matrix with values > 0.5

%%
%  plot 1D function
clear all;clc
t=0:0.01:10;
x=rand(1,length(t));
y1 = 0.5*cos(pi/2*t);
y2 = x.*cos(pi/2*t);

figure
pict=plot(t,y1,'--',t,y2,'-')
xlabel('time(s)')
ylabel('signal(V)')
legend('0.5*cos(x)','rand*cos(x)')
title('Example of multiple plots')
axis([0 2*pi -3 3])
%%
close all;clc
plot(t,y1)
hold on
scatter(t,y2)
%%
clear;clc
close all;  %2D PLOT

[X,Y] = meshgrid(-8:.5:8,-8:.5:8);
R = sqrt(X.^2 + Y.^2) + eps;
Z = sin(R)./R;
figure(1)
subplot(1,3,1)
contour(X,Y,Z)
subplot(1,3,2)
contourf(X,Y,Z)
subplot(1,3,3)
pcolor(X,Y,Z)


figure(2)
subplot(2,2,1)
plot3(X,Y,Z)
subplot(2,2,2)
mesh(X,Y,Z)
subplot(2,2,3)
surf(X,Y,Z)
subplot(2,2,4)
contour3(X,Y,Z)
%%
close all

%%
%function
clear all; clc
astring='today you are learning MATLAB';
number=123456;
[occur,locations]=searchit(astring,'a')
fprintf('number is still %d \n',number)



%%
% time your code
clear all;clc
tic
pause(3)
toc

%%
% if end
clear all; clc
flag=19;
if flag < 5
    disp('OK');
elseif flag < 9
    disp('Great');
else
    disp('Yes');
end
%%
clc
flag=9;
switch flag
    case  5
        disp('OK')
    case  9
        disp('Great')
end



%%
% for loop
clear all; clc;
for n=1:2:10
    disp(n)
end

%%
% while loop
clear all; clc %clean screen
n=1;
while n < 10
    n=n+1
end
n

%%
%break from a loop
clc %clean screan
for n=1:10
    n
    if n == 5
        break;
    end
end
n

%%
%continue in loop
clc  % clean screen
for n=1:10
    if n==5
        continue;
    end
    disp(n);
end

%%
%return in script/function
clc  % clean screen
val1=6;
val2=8;
multiply(val1,val2)


%%
% save and load data
clear all;clc
A=rand(10);
B=eye(10);
C=rand(10,1);
save myjob
clear all
who
load myjob
who
save myjob2 A
clear all
load myjob2
who
save myjob3.txt A -ascii
clear all;
%%
clc
text = fileread('seq.txt');  %load text from file
[count,array]=searchit(text,'GG')

%%
% stay away from for loop example
% given a matrix, find the sum of elements > 0.5

% for loop
clc
clear all;
L=10000;
a=rand(L);
tic
total=0;
for m=1:L
    for n=1:L
        if a(m,n) > 0.5
            total=total+a(m,n);
        end
    end
end
total
toc

%%
% vector version

tic
b=a>0.5;
c=a.*b;
total2=sum(sum(c))
toc



