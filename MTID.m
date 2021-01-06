close all
clear all

%单载频矩形脉冲信号
fs = 20e6; %sample rate
pri = 2.5e-5; %pulse width
t_single = 0:1/fs:(pri-1/fs);
ut_single = sin(20*pi/pri.*t_single);

%多个脉冲信号
fs_pulse = 20e6 %sample rate
T = 2.5e-4 %repeat period
n_pulse = 600;
t_pulse = 0:1/fs_pulse:(T-1/fs_pulse); %single pulse repeat time
t_npulse = 0:1/fs_pulse:(n_pulse*T-1/fs_pulse); %total time
utn_single1 = [ut_single,zeros(1,length(t_pulse)-length(ut_single))];
utn_single = repmat(utn_single1,1,n_pulse);
utn_single = utn_single(:)';

% 静止回波
 fs_decay = 20e6; %sample rate
 delta_t = 3e-5; %decay time
 n_delta = delta_t*fs_decay; %delay points
 A_decay = 1; %decay
 ut_dec = A_decay*ut_single;
 utn_1dec=[zeros(1,n_delta),ut_dec,zeros(1,length(t_pulse)-n_delta-...
length(ut_dec))];
utn_dec = repmat(utn_1dec,1,n_pulse);
utn_dec = utn_dec(:)';

% 运动回波 %fd = 0.5/T;
fd = 2000;
tf_shift = cos(2*pi*fd*t_npulse);
utn_target1=[zeros(1,n_delta),ut_dec,zeros(1,length(t_pulse)-n_delta-...
    length(ut_dec))];
utn_target = repmat(utn_target1,1,n_pulse);
utn_target = utn_target(:)';
utn_targetmov = tf_shift.*utn_target;

% 瑞利杂波叠加
load rayleigh_clutter_20dB.mat %加载杂波数据
ydata = ydata/40;
r_clutter2 = []; %瑞利杂波
for k = 1:floor(length(utn_targetmov)/(length(ydata)/2))
    start_c = floor(rand()*(length(ydata)/3))+1;
    r_clutter2 = [r_clutter2,ydata(start_c:start_c+(length(ydata)/2)-1)];
end
utn_cluttermov = utn_targetmov + r_clutter2;
figure
plot(t_npulse,utn_cluttermov);
title('瑞利杂波下单载频矩形脉冲运动回波信号');

% 叠加瑞利杂波后时频图
data_clumov = [];
for k = 0:n_pulse-1
    data_clumov=[data_clumov;utn_cluttermov(k*length(utn_1dec)+1:(k+1)*length(utn_1dec))];
end
n = 1:n_pulse;
figure
mesh(t_pulse*fs_pulse,n,abs(data_clumov));
colorbar
xlabel('time(s)');ylabel('frequency(Hz)');title('瑞利杂波下运动回波信号 STFT Time-Frequency Map');

% 运动目标 MTI 处理
n_accu2 = 4; %相消个数
mti_pulse=[1 -4 6 -4 1]; %相消器响应
data_accumov = [];
for i=1:length(utn_target1)
    temp=conv(mti_pulse,data_clumov(:,i)');
    data_accumov=[data_accumov temp(1,n_accu2/2:n_accu2/2-...
        1+length(data_clumov(:,i)))'];
end
figure
mesh(t_pulse*fs_pulse,n,abs(data_accumov));
colorbar;
xlabel('time/s');
ylabel('frequency/Hz');
title('MTI');

 % 运动目标 MTD 处理
 mtd_time_freq=[];
 for i=1:length(utn_target1)
     mtd_time_freq=[mtd_time_freq abs(fft(data_accumov(:,i)'))'];
 end
tmp=(0:round(1/pri)/length(data_accumov(:,i)):round(1/pri)-1);
figure
mesh(t_pulse*fs_pulse,tmp,abs(mtd_time_freq));
xlabel('距离单元');ylabel('多普勒频移');title('MTD');


 