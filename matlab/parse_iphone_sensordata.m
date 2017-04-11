function [ sensor_data ] = parse_iphone_sensordata( filename )
% Function to parse iPhone sensor data file from the app 'Data Collection'.
% Input: filename as a string
% Output: cell array with Attitude, Accelerometer, Gyro and Magnetometer data
%         in the form [timestamp, datapoint].

fid = fopen(filename,'r');
bob = textscan(fid,'%s %n %n %n %n','Headerlines',2);
fclose(fid);

npoints = length(bob{1})/4;
attitude_data = zeros(npoints,4);
accel_data = zeros(npoints,4);
gyro_data = zeros(npoints,4);
magneto_data = zeros(npoints,4);

attitude_inds = [1:4:(length(bob{1})-3)];    
accel_inds = [2:4:(length(bob{1})-2)];
gyro_inds = [3:4:(length(bob{1})-1)];
magneto_inds = [4:4:(length(bob{1}))];

%Attitude data
for i = 1:4

    attitude_data(:,i) = bob{i+1}(attitude_inds);

end
sensor_data{1} = attitude_data;


%Accelerometer data
for i = 1:4

    accel_data(:,i) = bob{i+1}(accel_inds);

end
sensor_data{2} = accel_data;


%Gyro data
for i = 1:4

    gyro_data(:,i) = bob{i+1}(gyro_inds);

end
sensor_data{3} = gyro_data;


%Magneto data
for i = 1:4

    magneto_data(:,i) = bob{i+1}(magneto_inds);

end
sensor_data{4} = magneto_data;

%% Plot accelerometer data
figure
plot(accel_data(:,1),accel_data(:,2),'b')
hold on
plot(accel_data(:,1),accel_data(:,3),'g')
plot(accel_data(:,1),accel_data(:,4),'r')
title('Accelerometer data')
xlabel('Timepoint')
ylabel('g')
legend('x','y','z')

%% Plot gyro data
figure
plot(gyro_data(:,1),gyro_data(:,2),'b')
hold on
plot(gyro_data(:,1),gyro_data(:,3),'g')
plot(gyro_data(:,1),gyro_data(:,4),'r')
title('Gyro data')
xlabel('Timepoint')
ylabel('rad s^{-1}')
legend('x','y','z')

%% Plot magneto data
figure
plot(magneto_data(:,1),magneto_data(:,2),'b')
hold on
plot(magneto_data(:,1),magneto_data(:,3),'g')
plot(magneto_data(:,1),magneto_data(:,4),'r')
title('Magnetometer data')
xlabel('Timepoint')
ylabel('\muT')
legend('x','y','z')

end

