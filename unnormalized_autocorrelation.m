clear
clc
close all

%for i = 1
    
%ch = int2str(i);
%data = strcat('Track(',ch,').txt');
%f_location = strcat('C:\Users\Bishwa Ranjan Si\Desktop\analysis_data\active+flow_data(matlab) _mod\0.5wt%\100 ul_hr\mov1(1.1.19)\',data);
%[n x y] = textread(f_location, '%f %f %f %*f %*f'); %read data from the text file

filename = '48m_80wt%_long_trajectory.xlsx';
n= xlsread(filename,'A:A');
x= xlsread(filename,'B:B');
y= xlsread(filename,'C:C');
x = x.*(1.66);
y = y.*(1.66);
t=zeros(length(n),1);

for i=1:length(n)
    t(i)=(n(i)-n(1))/10;
end

% x=smoothdata(x1(:),'rloess',6);
% y=smoothdata(y1(:),'rloess',6);
vx=zeros(length(x)-1,1);
vy=zeros(length(y)-1,1);
% limit=length(n)-1;
% A=0;
% klimit=round(limit/10,0);

%velocity calculation--------------------------------------
for i=1:length(x)-1
    vx(i)=(x(i+1)-x(i))/(t(i+1)-t(i));
    vy(i)=(y(i+1)-y(i))/(t(i+1)-t(i));
end
vt_old=[vx,vy];
%-----------------------------------------------------------

v_len=length(vx);
vt=zeros(v_len,2);
c=0;
j=1;
for i=1:1:v_len
    if (vt_old(i,1)~=0 && vt_old(i,2)~=0)
        vt(j,1)=vt_old(i,1);
        vt(j,2)=vt_old(i,2);
        j=j+1;
    else
        c=c+1;
    end
end

len_new=v_len-c;
vt=vt(1:len_new,:);

limit=len_new-1;
A=0;
klimit=round(limit/5,0);

%Autocorrelation calculation--------------------------------
for k=1:1:klimit
    for j=1:limit-k  
        A=A+(vt(j-1+k,1)*vt(j,1)+vt(j-1+k,2)*vt(j,2));
    end
 Ct(k,:)=A/(limit-k);
 A=0;
 delt(k,:)=k/10;
end

% plot(delt, Ct)
% plot(x,y)
% hold
% plot(x1,y1)
% hold off

% % %Taur Calculation---------------------------------------------
f = @(b,delt) b(1)*exp(-2*b(2).*delt)                                    % Objective Function
B = fminsearch(@(b) norm(Ct - f(b,delt)), [1 10])                 % Estimate Parameters
figure
plot(delt, Ct, 'pg')
% xlim([0 klimit/10])
% ylim([0 1])
hold on
plot(delt, f(B,delt), '-r')
hold off
grid
xlabel('delta t')
ylabel('C(t)')
