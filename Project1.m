%Clear all variable and close other files

clc;
clear all;
close all;

%User Defined Target Parameters  (Defined as per test case provided)

%Target 1

target1_distance = 300;                                     %Distance of the target from the radar in meters

target1_velocity_kmph = -50;                                %Velocity in kilometer per hour

target1_theta = 60;                                         %Angle of target in degrees (-ve angle indicates target left of the radar)

%Target 2

target2_distance = 200;                                     %Distance of the target from the radar in meters

target2_velocity_kmph = 75;                                 %Velocity in kilometer per hour

target2_theta = -45;                                        %Angle of target in degrees (-ve angle indicates target left of the radar)



%%

%Velocity unit conversion from km/hr to m/sec

target1_velocity = target1_velocity_kmph/3.6;               %Velocity of target 1 in m/sec 
target2_velocity = target2_velocity_kmph/3.6;               %Velocity of target 2 in m/sec 


%Constant Definitions

c = physconst('LightSpeed');                                %Speed of light in m/sec

%Radar Parameter Definitions

%Range parameters

range_max = 300;                                            %Maximum range of the radar in meters
range_res = 0.3;                                            %Desired range resolution in meters

%Velocity Parameters

vel_max = 100/3.6;                                          %Maximum velocity to be detected in m/sec (100 kmph) 
vel_res = 0.8/3.6;                                          %Desired velocity resolution in m/sec


%Angle Parameters

angle_res = 4.1;                                              %Desired minmum angle to be detected [From 1] for automotive applications

num_antennas = 2.^(nextpow2(floor((max(target1_distance, target2_distance))/(2*angle_res))+ 1)); %From [1]

%Chirp Signal Parameters (Computed)

total_signal_time = 0.03456;                                %Total observation time from [1]

num_chirps = 2.^(nextpow2(2*vel_max/vel_res));              %Number of chirps required for meeting desired velocity maximum and resolution (Taking powers of 2 for faster FFT computation)

num_samples_chirp = 2.^(nextpow2(range_max/range_res));     %Number of samples in each chirp (Taking powers of 2 for faster FFT computation)

chirp_time = total_signal_time/num_chirps;

chirp_BW = c/(2*range_res);                                 %From [4]

lambda_carrier = vel_res * 2 * num_chirps * chirp_time;     %From [4]

freq_carrier = c/lambda_carrier;                            %Carrier Frequency in GHz

k_chirp = chirp_BW/chirp_time;                              %Modulation index of chirps from [1]

f_beat_max = abs(2*vel_max/lambda_carrier) + k_chirp * 2 * range_max/c; %Maximum Beat Frequency (For finding Nquist frequency and thus sampling time) from [4]

sampling_frequency = 2.5 * f_beat_max;                      %Sampling Frequency based on Nyquist Theorem

sampling_time = 1/sampling_frequency;


%Beat Signal Synthesis

d1 = target1_distance;

v1 = target1_velocity;

d2 = target2_distance;

v2 = target2_velocity;

ds = lambda_carrier/2;                                      %Distance between Antennas [1] mentions Dan as lamda-carrier however that limits the maximum angle to 30 degrees

BW = chirp_BW;

Ts = sampling_time;

Tchirp = chirp_time;

n = linspace(0,1,num_samples_chirp);                        

beat_signal = zeros(num_chirps, num_samples_chirp,num_antennas);        %Preallocation of size of 2D matrix for data


%Filling the 3D data matrix with chirp_signal x samples of chirp signal

for k = 1: num_antennas

    for l=1:num_chirps

        for n1=1:length(n)

            beat_signal(l,n1,k) = exp(1j*2*pi*((2*BW*d1*n1*Ts)/(Tchirp*c) + 2*v1*(l-1)*Tchirp/lambda_carrier + ds*sind(target1_theta)*(k-1)/lambda_carrier)) + exp(1j*2*pi*((2*BW*d2*n1*Ts)/(Tchirp*c) + 2*v2*(l-1)*Tchirp/lambda_carrier + ds*sind(target2_theta)*(k-1)/lambda_carrier)); %Beat signal synthesis using exponential
        end
    end

end

%2D FFT for Range-Doppler Velocity 

beat_signal_fft = fft2(beat_signal(:,:,1), num_chirps, num_samples_chirp);          %FFT computed for first 2D Beat Signal Matrix

beat_signal_fft_mag = abs(beat_signal_fft);

beat_signal_fft_mag_shifted = abs(fftshift(beat_signal_fft));                       %FFT Shift to make the FFT zero centered


%Axis Scaling for frequency bins to correct range and velocity units

x_axis_unscaled = linspace(-0.5,0.5,num_samples_chirp);

x_axis = (x_axis_unscaled)*(chirp_time*c)/(2*BW*sampling_time);                     %From [4]

y_axis_unscaled = linspace(-0.5,0.5, num_chirps);                                   %From [4]

y_axis = y_axis_unscaled * 3.6 * (lambda_carrier)/(2*chirp_time);                   %From [4]

doppler_range_image = imagesc(x_axis, y_axis, beat_signal_fft_mag_shifted); 

set(gca,'YDir','normal')                                                            %Flips the y-axis! (to see increasing values, not decreasing)
colorbar                                                                            %Display the Colour bar   
colormap jet
title("Distance vs Doppler Velocity");
xlabel("Distance (in Meters)")
ylabel("Doppler Velocity (in km/hr)")


%Finding the max values in the Distance - Doppler Velocity FFT 2D matrix

max_value1 = 0;                                                                     %Initialization of first maximum value

%For loop for finding first maximum value

for row = 1:num_chirps
    for col = 1:num_samples_chirp

        if beat_signal_fft_mag_shifted(row,col) > max_value1

            max_value1 = beat_signal_fft_mag_shifted(row,col);
        end
    end

end

[max1_row_index, max1_column_index] = find(beat_signal_fft_mag_shifted == max_value1);  %Index for first maximum


%Conditions so that the second maximum value is not chosen from the same
%cluster thus resulting in missing one of the Targets

row_top = max1_row_index - 2;
row_bottom = max1_row_index + 2;

col_left = max1_column_index - 5;
col_right = max1_column_index + 5;

%For loop for finding second maximum value (To add the condition that peaks are not from the same cluster)

max_value2 = 0;

for row=1:num_chirps
    for col=1:num_samples_chirp
        if (((1 <= row_top) && (row_top <= num_chirps))||((1 <= row_bottom))  && ((row_bottom <= num_chirps))) && (((1 <= col_left) && (col_left <= num_samples_chirp))||((1 <= col_right) && (col_right <= num_samples_chirp)))
        
            if (row < row_top) || (row > row_bottom)

                if (col < col_left) || (col > col_right)

                    if (beat_signal_fft_mag_shifted(row,col) > max_value2) &&  (beat_signal_fft_mag_shifted(row,col) ~= max_value1)

                        max_value2 = beat_signal_fft_mag_shifted(row,col);

                    end
                end
            end
        end
    end
end

[max2_row_index, max2_column_index] = find(beat_signal_fft_mag_shifted == max_value2);              %Index for second maximum

%Plotting 1D FFT Slice for Range along the row with the max amplitude

x_axis_unscaled_1d = linspace(-0.5, 0.5,num_samples_chirp);

Xaxis_distance_1d = (x_axis_unscaled_1d)*(chirp_time*c)/(2*BW*sampling_time);

normalized_distance_1 = beat_signal_fft_mag_shifted(max1_row_index,:)/max(beat_signal_fft_mag_shifted(max1_row_index,:));  %FFT magnitude is normalized to view the peaks more accurately

normalized_distance_2 = beat_signal_fft_mag_shifted(max2_row_index,:)/max(beat_signal_fft_mag_shifted(max2_row_index,:));

%Plot for Distance using 1D FFT Slice

figure

plot(Xaxis_distance_1d, normalized_distance_1,'r', Xaxis_distance_1d, normalized_distance_2,'b');  

title("Distance Plot");
xlabel("Distance (in meters)");
ylabel("Normalized FFT Magnitude");

%Distance Detected

distance_targets_estimated = [Xaxis_distance_1d(max1_column_index), Xaxis_distance_1d(max2_column_index)];


%Plotting 1D FFT Slice for velocity along the column with the max amplitude

Xaxis_velocity_unscaled_1d = linspace(-0.5, 0.5,num_chirps);

Xaxis_velocity_1d = Xaxis_velocity_unscaled_1d * 3.6 * (lambda_carrier)/(2*chirp_time);

normalized_velocity_1 = beat_signal_fft_mag_shifted(:,max1_column_index)/max(beat_signal_fft_mag_shifted(:,max1_column_index)); %FFT magnitude is normalized to view the peaks more accurately

normalized_velocity_2 = beat_signal_fft_mag_shifted(:,max2_column_index)/max(beat_signal_fft_mag_shifted(:,max2_column_index));

%Plot for Velocity using 1D FFT Slice

figure

plot(Xaxis_velocity_1d,normalized_velocity_1,'r', Xaxis_velocity_1d, normalized_velocity_2,'b');  

title("Doppler Velocity Plot");
xlabel("Velocity in (km/hr)");
ylabel("Normalized FFT Magnitude");

%Velocity Detected

velocity_targets_estimated = [Xaxis_velocity_1d(max1_row_index), Xaxis_velocity_1d(max2_row_index)];


%For finding the Azimuth (Angle of arrival for the targets)

%Computing FFT with the 3rd Dimension 

azimuth_fft = fftn(beat_signal(:,:,:), [num_chirps, num_samples_chirp, num_antennas*32]);     %Zero padding is performed for a better looking FFT plot

azimuth_fft_mag = abs(fftshift(azimuth_fft));

azimuth_fft_target1 = azimuth_fft_mag(max1_row_index, max1_column_index, :);                  %Max value 1 indices are taken from Distance-Velocity 2D FFT for ploting azimuth information

azimuth_fft_target2 = azimuth_fft_mag(max2_row_index, max2_column_index, :);                  %Max value 1 indices are taken from Distance-Velocity 2D FFT for ploting azimuth information

normalized_azimuth_fft_target1 = azimuth_fft_mag(max1_row_index, max1_column_index, :)/max(azimuth_fft_mag(max1_row_index, max1_column_index, :));  %Normalized FFT Magnitude to better view the peaks

normalized_azimuth_fft_target2 = azimuth_fft_mag(max2_row_index, max2_column_index, :)/max(azimuth_fft_mag(max2_row_index, max2_column_index, :)); 

Xaxis_angle_unscaled = asind(linspace(-1,1,length(azimuth_fft_target1)));

%Plot for Azimuth Angle

figure

plot(Xaxis_angle_unscaled ,reshape(normalized_azimuth_fft_target1,[num_antennas*32,1]), 'r', Xaxis_angle_unscaled, reshape(normalized_azimuth_fft_target2, [num_antennas*32,1]), 'b');

title("Azimuth","FontSize",15);
xlabel("Azimuth Angle (in Degrees)");
ylabel("Normalized FFT Magnitude");


%For Azimuth-Range 2D Doppler plot (If velocity information is not required)

beat_signal_1 = beat_signal(1,:,:);                                         %Considering the first chirp ignal across all receiver antennas

beat_signal_1 = reshape(beat_signal_1, [num_samples_chirp, num_antennas]);

beat_signal_1_fft = fft2(beat_signal_1, num_samples_chirp*8, num_antennas*32); %Zero padding for better viewing resolution

beat_signal_1_fft_mag = abs(fftshift(beat_signal_1_fft));

x_axis_azi = linspace(-1,1, num_samples_chirp*8);

y_axis_azi = (x_axis_unscaled)*(chirp_time*c)/(2*BW*sampling_time);

%2D Plot for Distance and Azimuth Angle information

figure

azimuth_range_image = imagesc(x_axis_azi, y_axis_azi, beat_signal_1_fft_mag); 

set(gca,'YDir','normal')                                                            %Flips the y-axis! (to see increasing values, not decreasing)
colorbar                                                                            %Display the Colour bar   
colormap jet
title("Distance vs Sine(Azimuth)");
xlabel("Sine(Azimuth) (in Radians)")
ylabel("Distance (in Meters)")

%Showing Range and Azimuth using 1D FFT Slices

%Finding the maximum values along the rows for the distance

max_value_azi1 = 0; %First maximum value

%For loop for finding first maximum value for Azi-Range

for row_azi = 1:num_samples_chirp*8
    for col_azi = 1:num_antennas*32

        if beat_signal_1_fft_mag(row_azi,col_azi) > max_value_azi1

            max_value_azi1 = beat_signal_1_fft_mag(row_azi,col_azi);

        end
    end

end

[max1_azi_row_index, max1_azi_column_index] = find(beat_signal_1_fft_mag == max_value_azi1); %Index of max value 1

%Condition for max value 2 so that the same targets are not selected

row_top_azi = max1_azi_row_index - 10;
row_bottom_azi = max1_azi_row_index + 10;

col_left_azi = max1_azi_column_index - 30;
col_right_azi = max1_azi_column_index + 30;

%For loop for finding second maximum value

max_value_azi2 = 0;

for row_azi=1:num_samples_chirp*8
    for col_azi=1:num_antennas*32
        if (((1 <= row_top_azi) && (row_top_azi <= num_samples_chirp*8))||((1 <= row_bottom_azi))  && ((row_bottom_azi <= num_samples_chirp*8))) && (((1 <= col_left_azi) && (col_left_azi <= num_antennas*32))||((1 <= col_right_azi) && (col_right_azi <= num_antennas*32)))
        
            if (row_azi < row_top_azi) || (row_azi > row_bottom_azi)

                if (col_azi < col_left_azi) || (col_azi > col_right_azi)

                    if (beat_signal_1_fft_mag(row_azi,col_azi) > max_value_azi2) &&  (beat_signal_1_fft_mag(row_azi,col_azi) ~= max_value_azi1)

                        max_value_azi2 = beat_signal_1_fft_mag(row_azi,col_azi);

                    end
                end
            end
        end
    end
end

[max2_azi_row_index, max2_azi_column_index] = find(beat_signal_1_fft_mag == max_value_azi2);


%Plotting 1D FFT Slice for Azimuth Range along the column with the max amplitude

x_axis_unscaled_1d_azi = linspace(-0.5, 0.5,num_samples_chirp*8);

Xaxis_distance_1d_azi = (x_axis_unscaled_1d_azi)*(chirp_time*c)/(2*BW*sampling_time);

normalized_distance_1_azi = beat_signal_1_fft_mag(:,max1_azi_column_index)/max(beat_signal_1_fft_mag(:,max1_azi_column_index));

normalized_distance_2_azi = beat_signal_1_fft_mag(:,max2_azi_column_index)/max(beat_signal_1_fft_mag(:,max2_azi_column_index));

%Plot for Range using Azimuth-Range Dimension

figure

plot(Xaxis_distance_1d_azi, normalized_distance_1_azi,'r', Xaxis_distance_1d_azi, normalized_distance_2_azi,'b');  

title("Distance Plot using Aziumth-Distance Dimension");
xlabel("Distance (in meters)");
ylabel("Normalized FFT Magnitude");


%Plotting 1D FFT Slice for Azimuth Angle along the row with the max amplitude

x_axis_1d_azi_angle = asind(linspace(-1, 1,num_antennas*32)); %In degrees

normalized_angle_1_azi = beat_signal_1_fft_mag(max1_azi_row_index,:)/max(beat_signal_1_fft_mag(max1_azi_row_index,:));

normalized_angle_2_azi = beat_signal_1_fft_mag(max2_azi_row_index,:)/max(beat_signal_1_fft_mag(max2_azi_row_index,:));

%Plot for Azimuth angle using Azimuth-Range Dimension

figure

plot(x_axis_1d_azi_angle, normalized_angle_1_azi,'r', x_axis_1d_azi_angle, normalized_angle_2_azi,'b');  

title("Azimuth Angle Plot using Aziumth Dimension");
xlabel("Angle (in Degrees)");
ylabel("Normalized FFT Magnitude");


%Azimuth Angle Detected

azi_angle_estimated = [x_axis_1d_azi_angle(max1_azi_column_index), x_axis_1d_azi_angle(max2_azi_column_index)];


%Birds Eye Plot

mag_x1 = distance_targets_estimated(1)*sind(azi_angle_estimated(2));

mag_y1 = distance_targets_estimated(1)*cosd(azi_angle_estimated(2));

mag_x2 = distance_targets_estimated(2)*sind(azi_angle_estimated(1));

mag_y2 = distance_targets_estimated(2)*cosd(azi_angle_estimated(1));

figure

plot(mag_x1, mag_y1,mag_x2, mag_y2,'Marker', '*');

axis([-range_max, range_max, 0, range_max]);

grid on
title("Birds Eye view of the targets")
xlabel("Lateral Distance (in Meters)")
ylabel("Distance in Front")
legend('Radar at 0,0')



%Displaying the estimated distance, velocity and angle information:

message_line1 = "Detected Target Information"; 

disp(message_line1);


fprintf(1, '\n');

message_line2 = ['Distance: ', num2str(distance_targets_estimated(2)), ' meters   and ', num2str(distance_targets_estimated(1)),' meters'];

message_line3 = ['Velocity: ', num2str(velocity_targets_estimated(2)),' km/hr   and ', num2str(velocity_targets_estimated(1)),' km/hr'];

message_line4 = ['Azimuth Angle: ', num2str(azi_angle_estimated(1)),' degrees   and ', num2str(azi_angle_estimated(2)),' degrees'];

disp(message_line2);

fprintf(1, '\n');

disp(message_line3);


fprintf(1, '\n');

disp(message_line4);


fprintf(1, '\n');
