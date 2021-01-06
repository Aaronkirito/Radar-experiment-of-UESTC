clear all
close all
% LFM�����ź�
LFM_B = 2e8; %200MHz
RE_t = 1e-6; %������
LFM_fs = 4*LFM_B; %LFM�źŲ���Ƶ��
LFM_A = 10; %LFM�źŷ���
LFM_f0 = 1000; %�ز�Ƶ��
point = floor(RE_t/(1/LFM_fs));%��������
LFM_N = linspace(-RE_t/2,RE_t/2,point); %������
y_LFM = LFM_A * exp(1j*(2*pi*LFM_f0*LFM_N+pi*(LFM_B/RE_t)*LFM_N.^2));
figure('name','LFM�ź�');
subplot(2,1,1)
plot(1e6*LFM_N,y_LFM);
xlabel('t/us'); title('LFM�źŲ���');
subplot(2,1,2)
freq = linspace(-LFM_fs/2,LFM_fs/2,point);%Ƶ�������
sf = fftshift(fft(y_LFM));
plot(freq/1e6,abs(sf));
xlabel('f/MHz');title('�ز��źŻ���Ƶ��');

%�Ϳ����ź�
BPSK_fs = 20e6; %������
B_code = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 -1];
num = length(B_code); %�Ϳ���λ��
BPSK_B = 13e-6; %����
BPSK_tao = 0:1/BPSK_fs:(BPSK_B-1/BPSK_fs); %��λ����ʱ��
BPSK_A = repelem(B_code,floor(length(BPSK_tao)/num)); %��ֵ
BPSK_y = BPSK_A.*sin(2*pi/(BPSK_B/num).*BPSK_tao); %BPSK�����ź�
BPSK_N = linspace(0,BPSK_fs,length(BPSK_y));
BPSK_fft = abs(fft(BPSK_y));
figure('name','BPSK�ź�')
subplot(2,1,1)
plot(BPSK_tao,BPSK_y);
title('BPSK�źŲ���');
subplot(2,1,2)
plot(BPSK_N,BPSK_fft);
title('�ز��źŻ���Ƶ��');

%��Ƶ�ź�
dan_fs = 300; %������
dan_num = 600;
dan_seq = (0:dan_num-1)*(1/dan_fs); %ʱ������
dan_A = 10;
dan_f = 80;
dan_y = dan_A*cos(2*pi*dan_f*dan_seq);
figure('name','��Ƶ�ź�');
subplot(2,1,1)
plot(dan_seq,dan_y);
xlabel('t/s'); ylabel('����/V'); title('��Ƶ�źŲ���');
dan_N = 2^nextpow2(dan_num); %��������
dan_A = fft(dan_y,dan_N)/dan_N*2 %��ֵ
dan_f = dan_fs/dan_N*(0:1:dan_N-1); %Ƶ��
dan_amp = abs(dan_A);
subplot(2,1,2)
plot(dan_f(1:dan_N/2),dan_amp(1:dan_N/2));
xlabel('Ƶ��/Hz'); ylabel('��ֵ');title('�ز��źŻ���Ƶ��');

%��������ź�
fs_pulse = 20e6; %������
RE_t = 1e-4; %�ظ�����
n_pulse = 10; %�������
n_delta = 600; %�ӳٸ���
A_decay = 0.6; %˥��
v_target = 2000;
lambda = 3e8 * 1e-8; %����
fd = 2*v_target/lambda; %������Ƶ��
t_pulse = 0:1/fs_pulse:(RE_t-1/fs_pulse); %�������ظ�ʱ��
t_npulse = 0:1/fs_pulse:(n_pulse*RE_t-1/fs_pulse); %��n������ʱ��
tf_shift = cos(2*pi*fd*t_npulse);

%LFM�źŻز�ģ��
LFM_utn_1 = [y_LFM,zeros(1,length(t_pulse)-length(y_LFM))];
LFM_utn = repmat(LFM_utn_1,1,n_pulse);
LFM_utn = LFM_utn(:)';
figure('name','LFM �ز�ģ��')
subplot(3,1,1)
plot(t_npulse,LFM_utn)
title('LFM �ź�');
LFM_ut_dec = A_decay*y_LFM;
LFM_utn_1dec=[zeros(1,n_delta),LFM_ut_dec,zeros(1,length(t_pulse)-n_delta-length(LFM_ut_dec))];
LFM_utn_dec = repmat(LFM_utn_1dec,1,n_pulse);
LFM_utn_dec = LFM_utn_dec(:)';
subplot(3,1,2)
hold on
plot(t_npulse,LFM_utn);
plot(t_npulse,LFM_utn_dec);
hold off
title('LFM ��ֹ�ز��ź�');
LFM_utn_target1 = [zeros(1,n_delta),LFM_ut_dec,zeros(1,length(t_pulse)-n_delta-length(LFM_ut_dec))];
LFM_utn_target = repmat(LFM_utn_target1,1,n_pulse);
LFM_utn_target = LFM_utn_target(:)';
LFM_utn_targetmov = tf_shift.*LFM_utn_target;
subplot(3,1,3)
hold on
plot(t_npulse,LFM_utn);
plot(t_npulse,LFM_utn_targetmov);
hold off
title('LFM�˶��ز��ź�');

%BPSK�źŻز�ģ��
BPSK_utn_1 = [BPSK_y,zeros(1,length(t_pulse)-length(BPSK_y))];
BPSK_utn = repmat(BPSK_utn_1,1,n_pulse);
BPSK_utn = BPSK_utn(:)';
figure('name','BPSK �źŻز�ģ��')
subplot(3,1,1)
plot(t_npulse,BPSK_utn);
title('BPSK �ź�');
BPSK_ut_dec = A_decay*BPSK_y;
BPSK_utn_1dec = [BPSK_ut_dec,zeros(1,length(t_pulse)-length(BPSK_ut_dec))];
BPSK_utn_dec = repmat(BPSK_utn_1dec,1,n_pulse);
BPSK_utn_dec = BPSK_utn_dec(:)';
subplot(3,1,2)
hold on
plot(t_npulse,BPSK_utn);
plot(t_npulse,BPSK_utn_dec);
hold off
title('BPSK��ֹ�ز��ź�');
BPSK_utn_target1=[zeros(1,n_delta),BPSK_ut_dec,zeros(1,length(t_pulse)-n_delta-length(BPSK_ut_dec))];
BPSK_utn_target = repmat(BPSK_utn_target1,1,n_pulse);
BPSK_utn_target = BPSK_utn_target(:)';
BPSK_utn_targetmov = tf_shift.*BPSK_utn_target;
subplot(3,1,3)
hold on
plot(t_npulse,BPSK_utn);
plot(t_npulse,BPSK_utn_targetmov);
hold off
title('BPSK �˶��ز��ź�');

%��Ƶ�źŻز�ģ��
dan_utn_1 = [dan_y,zeros(1,length(t_pulse)-length(dan_y))];
dan_utn = repmat(dan_utn_1,1,n_pulse);
dan_utn = dan_utn(:)';
figure('name','��Ƶ�źŻز�ģ��')
subplot(3,1,1)
plot(t_npulse,dan_utn);
title('��Ƶ�ź�');
dan_ut_dec = A_decay*dan_y;
dan_utn_1dec=[zeros(1,n_delta),dan_ut_dec,zeros(1,length(t_pulse)-n_delta-length(dan_ut_dec))];
dan_utn_dec = repmat(dan_utn_1dec,1,n_pulse);
dan_utn_dec = dan_utn_dec(:)';
subplot(3,1,2)
hold on
plot(t_npulse,dan_utn);
plot(t_npulse,dan_utn_dec);
hold off
title('��Ƶ�źž�ֹ�ز��ź�');
dan_utn_target1=[zeros(1,n_delta),dan_ut_dec,zeros(1,length(t_pulse)-n_delta-length(dan_ut_dec))];
dan_utn_target = repmat(dan_utn_target1,1,n_pulse);
dan_utn_target = dan_utn_target(:)';
dan_utn_targetmov = tf_shift.*dan_utn_target;
subplot(3,1,3)
hold on
plot(t_npulse,dan_utn);
plot(t_npulse,dan_utn_targetmov);
hold off
title('��Ƶ�ź��˶��ز��ź�');

%LFM�ź�ģ��ͼ
LFM_tp = 1e-6;
LFM_B = 2e6;
point = 90;
[LFM_N, fd] = meshgrid(linspace(-1, 1, point), linspace(-20, 20, 2*point));
LFM_A = (1-abs(LFM_N*LFM_tp)/LFM_tp).*abs(sin(pi*(fd/LFM_tp*LFM_tp+LFM_B*LFM_N*LFM_tp).*(1-abs(LFM_N*LFM_tp)/LFM_tp))...
    /pi./(fd/LFM_tp*LFM_tp+LFM_B*LFM_N*LFM_tp)./(1-abs(LFM_N*LFM_tp)/LFM_tp));
A_norm = LFM_A/max(max(LFM_A));
figure('name','LFM �ź�ģ������ͼ');
mesh(LFM_N, fd, A_norm);
colormap jet;
title('LFM �ź�ģ������ͼ');
xlabel('ʱ��');
ylabel('������Ƶ��');

%BPSK�ź�ģ��ͼ
BPSK_N = 13; %2, 3, 4, 5, 7, 11, 13
    switch BPSK_N
        case 2 %�԰��ƽ˥��6dB
            if 1
                Barker_code = [1 -1];
            else
                Barker_code = [1 1];
            end
        case 3 %�԰��ƽ˥��9.5dB
            Barker_code = [1 1 -1];
        case 4 %�԰��ƽ˥��12.0dB
            if 1
                Barker_code = [1 1 -1 1];
            else
                Barker_code = [1 1 1 -1];
            end
        case 5 %�԰��ƽ˥��14.0dB
            Barker_code = [1 1 1 -1 1];
        case 7 %�԰��ƽ˥��16.9dB
            Barker_code = [1 1 1 -1 -1 1 -1];
        case 11 %�԰����˥��20.8dB
            Barker_code = [1 1 1 -1 -1 -1 1 -1 -1 1 -1];
        case 13 %�԰��ƽ˥��22.3dB
            Barker_code = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
        otherwise
            disp '�Ϳ��벻����';
            return;
    end
    
BPSK_tau = BPSK_N;
BPSK_samp_num = BPSK_N*10;
nfft = 2^ceil(log(BPSK_samp_num)/log(2));
u(1:nfft) = 0;
j = 0;
    for index = 1:10:BPSK_samp_num
        j = j+1;
        u(index:index+10-1) = Barker_code(j);
    end
    
delay = linspace(-BPSK_tau, BPSK_tau, nfft);
freq_del = 12/BPSK_tau/100;
j = 0;
vfft = fft(u,nfft);
    for freq = -6/BPSK_tau:freq_del:6/BPSK_tau
        j = j+1;
        exf = exp(1i*2*pi*freq*delay);
        u_times_exf = u.*exf;
        ufft = fft(u_times_exf, nfft);
        prod = ufft.*conj(vfft);
        ambig(:, j) = fftshift(abs(ifft(prod))');
    end
freq = -6/BPSK_tau:freq_del:6/BPSK_tau;
delay = linspace(-BPSK_N, BPSK_N, nfft);
value = 10*BPSK_N;
figure('name','�Ϳ���ģ������ͼ');
mesh(freq, delay, ambig/max(max(ambig)));
colormap jet;axis tight;
title('�Ϳ���ģ������ͼ');
xlabel('ʱ��');
ylabel('������Ƶ��');

%��Ƶ�ź�ģ��ͼ
dan_tp = 100e-6;
[t, fd] = meshgrid(linspace(-1, 1, point), linspace(-10, 10, 2*point));
dan_A = abs(sin(pi*fd/dan_tp*dan_tp.*(1-...
    abs(t*dan_tp)/dan_tp))/pi./fd/dan_tp/dan_tp./(1-abs(t*dan_tp)/dan_tp))...
    .*(1-abs(t*dan_tp)/dan_tp);
A_norm = dan_A/max(max(dan_A));
figure('name','��Ƶ�ź�ģ������ͼ');
mesh(t, fd, A_norm);
colormap jet;
title('��Ƶ�ź�ģ������ͼ');
xlabel('ʱ��');
ylabel('������Ƶ��');