function [b xdata ydata step msd_x msd_y] = MSD_old(data) %[tau,msd_x,msd_y] = MSD(data)
% INPUT
% x   = x coordinates of the particle trajectory in time
% y   = y coordinates of the particle trajectory in time
% t   = time, should be an integer value. Unit = 1/frames_per_second
% T   = Maximum time upto which mean squared displacement is to be
%       calculated, should be an integer value. Unit = 1/frames_per_second.
%       Greater T will result in less averaging for evaluating MSD(T)
% data = load('P1.txt');

filename = '19microns_only water.xlsx';
t = xlsread(filename,'A:A');
x = xlsread(filename,'B:B');
y = xlsread(filename,'C:C');

n = length(x);    % n = number of position data points
T = 2750;
% if T == 0         % T = 0, set the default value to maximum time possible
%     T = n - 1;
% end

l     = T;
msd_x = zeros(1,l);
std_x = zeros(1,l);
msd_y = zeros(1,l);
std_y = zeros(1,l);
tau   = zeros(1,l); % tau = time array corresponding to msd_x, msd_y arrays


for i = 1:1:l
    count  = 0;
    j      = 1;
    tau(i) = i;
   
    while j+i <= n % Exit if the end point of a time interval exceeds n
        % Locating next time origin if trajectory is discontinuous
        while j + i <= n && t(j + i) - t(j) > 1.1 * tau(i)
            % A factor of 1.1 is chosen because for the smallest
            % discontunity in the trajectory (i.e. one missing snapshot)
            % the time gap will be 2 units.
            j = j + 1;
        end
        % Evaluate MSD if the end point of the time interval has not
        % exceeded n
        if j + i <= n
            msd_x(i) = msd_x(i) + (x(j) - x(j + i))^2;
            msd_y(i) = msd_y(i) + (y(j) - y(j + i))^2;
            count    = count + 1;
            
            %[j,j+i,count,x(j),x(j+i),(x(j)-x(j+i))^2,y(j),y(j+i),(y(j)-y(j+i))^2]
            %pause
            j        = j + i;
        end
    end
    if count ~= 0
        msd_x(i) = msd_x(i)/count;
        msd_y(i) = msd_y(i)/count;
    end
end

%Calculate Standard deviation
for i = 1:1:l
    count  = 0;
    j      = 1;
    tau(i) = i;
   
    while j+i <= n % Exit if the end point of a time interval exceeds n
        % Locating next time origin if trajectory is discontinuous
        while j + i <= n && t(j + i) - t(j) > 1.1 * tau(i)
            % A factor of 1.1 is chosen because for the smallest
            % discontunity in the trajectory (i.e. one missing snapshot)
            % the time gap will be 2 units.
            j = j + 1;
        end
        % Evaluate MSD if the end point of the time interval has not
        % exceeded n
        if j + i <= n
            std_x(i) = std_x(i) + ((x(j) - x(j + i))^2 - msd_x(i))^2;
            std_y(i) = std_y(i) + ((y(j) - y(j + i))^2 - msd_y(i))^2;
            count    = count + 1;
            %[j,j+i,count,x(j),x(j+i),(x(j)-x(j+i))^2,y(j),y(j+i),(y(j)-y(j+i))^2]
            %pause
            j        = j + i;
        end
    end
    if count ~= 0
        std_x(i) = sqrt(std_x(i)/(count - 1));
        std_y(i) = sqrt(std_y(i)/(count - 1)); 
    end
end

%Change units from pixels^2 to um^2

msd_x = msd_x *(1)^2;
msd_y = msd_y *(1)^2;


disp ('diffusivities calculated from single point');
b(1,1) = msd_x(1)/0.208333/2;
b(1,2) = msd_y(1)/0.208333/ 2;
%b(2,1) = 0.5 * sum(tau(:).*msd_x(:))/sum(tau(:).*tau(:))*0.187^2 * 30;
%b(2,2) = 0.5 * sum(tau(:).*msd_y(:))/sum(tau(:).*tau(:))*0.187^2 * 30;

xdata(:,1) = msd_x;
ydata(:,1) = msd_y;
step(:,1) = tau*(1/29);

fprintf('Diff_X = %6.3g; Diff_Y = %6.3g\n', b(1,1), b(1,2));

plot(step,msd_x,'ob',step,msd_y,'sr');
xlabel('\tau');
ylabel('Mean Squared Displacement');
legend('msd_x','msd_y');
%saveas(gca, data(1:end-4),'jpg');

%linear fitting and plotting
px=polyfitZero(transpose(step),msd_x,1);
py=polyfitZero(transpose(step),msd_y,1);
pxfit = px(1)*step+px(2);
pyfit = py(1)*step+py(2);
hold on;
plot(step,pxfit,'b',step,pyfit,'r');

%diffusivities from linear regression
disp ('diffusivities calculated from linear fits');
Dx=px/2;
Dy=py/2;
fprintf('Dx = %6.3g; Dy = %6.3g\n', Dx(1), Dy(1));
% fidx = fopen(['MSDx-' data], 'a');
% fprintf(fidx, '%s\t', data);
% fprintf(fidx, '%f\t', msd_x');

filename = 'output.xlsx';
Sheet=strcat(['MSD-' data], 'a','data');
size=length(step);
A=zeros(size+2,3);
A = [step msd_x' msd_y'];
A(size+2,2)=Dx(1);
A(size+2,3)=Dy(1);
xlswrite(filename,A,Sheet);

%fprintf(fidx, '\n');
%fprintf(fidx, '%s\t', data);
%fprintf(fidx, '%f\t', std_x);
% fprintf(fidx, '\n');
% fclose(fidx);

% fidy = fopen(['MSDy-' data], 'a');
% fprintf(fidy, '%s\t', data);
% fprintf(fidy, '%f\t', msd_y');
% %fprintf(fidy, '\n');
% %fprintf(fidy, '%s\t', data);
% %fprintf(fidy, '%f\t', std_y);
% fprintf(fidy, '\n');
% fclose(fidy);
