
%========================================%
% 调频测距循环
%========================================%
clc
clear
close all

c = 3e8;
fs = 100e6;

% fr_c = 883500
 fr_c = 891200;
%--------全相位滤波预处理----------
Ns = 8192;
win = hanning(Ns); win1 = hann(Ns);
win2 = conv(win, win1);
win2 = win2 / sum(win2);

win_1 = win/sum(win);
winn = conv(win, win);
win_2 = winn/sum(winn);

filename = 'Data_34m_1.dat';  
fid = fopen(filename,'r');
% data_q = fread(fid,inf,'uint16');
N = 40000;
%--------以下是抗混叠滤波器-------
b = [0.0432446018265606, -0.0146985229518258, 0.157558806755585, -0.00627836872663652, 0.236193764779514, 0.0269245428655922, 0.220307816321058, 0.0124593189362773, 0.136131042329674, -0.00576003661136978, 0.0376979316762628;];
a = [1, -3.49404028503486, 8.63279997018137, -14.0423560473518, 17.9440410832339, -17.3105470156172, 13.2731977637244, -7.68419838946513, 3.34462242571084, -0.967053332362388, 0.155218614901453;];
a_BP= a.*exp(j*2*pi*20.5e6/fs*(0:10));
b_BP= b.*exp(j*2*pi*20.5e6/fs*(0:10));


for i = 1:1000
    fseek(fid,46592*2*2*(i-1),'bof');                   
    data = fread(fid,46594*2,'uint16');           
    
    data1 = data(1:2:end);   %参考波
    data2 = data(2:2:end);   %回波
    

    %figure;plot(1:45000,data1(1:45000),'r',1:45000,data2(1:45000),'b');
    m1 = data2(2001:42000);    
    r1 = data1(2001:42000);
    m1 = m1 - 8192;    %采样的原因是8000左右的，去掉均值
    r1 = r1 - 8192;
    
    x_r = filter(b_BP,a_BP,r1);    
    
    m_c = m1(1001:1000+Ns);
    f_m_c = fft(m_c,Ns);   %fft
    f_m_c(1:round(Ns/20)) = 0;
    f_m_c(Ns/2:end) = 0;
    [maxi_abs(i),m_posi(i)] = max(abs(f_m_c));
    %figure; plot(abs(f_m_c(0:2*m_posi(i))))
    fm_c = (m_posi(i)-1)/Ns*fs;
    Ratio = fm_c./fr_c;    
    
    angle_1 = angle(x_r);
    %figure; plot(angle_1)
    angle_2 = unwrap(angle_1);
    %figure; plot(angle_2)
    t = [0:N-1]*1./fs;
    angle1 = angle_2 - 2*pi*fr_c*t.';
    %figure;plot(angle1)
    
    
    a_SF= a_c.*exp(j*2*pi*fm_c/fs*(0:10));
    b_SF= b_c.*exp(j*2*pi*fm_c/fs*(0:10));
    
    x_m_c = filter(b_SF,a_SF,m1);
    %x_m_c(1:250) = 0;
    
    x_m_c = x_m_c.*exp(-j*Ratio*angle1);
    
     [maxi,m_posi_new(i)] = max(abs(fft(x_m_c)));
     %figure; plot(abs(fft(x_m_c)))
     fm_new1 = (m_posi_new(i)-1)./N*fs;
     angle1_s = angle1;
     x_m_c = x_m_c.*exp(-j*(fm_new1-fm_c)./fr_c*angle1_s);     
    

%-------俩个序列都做apfft
s1 = x_m_c(1 : 2*Ns-1);
    y1 = s1 .* win2;
    y1a = y1(Ns:end) + [0; y1(1:Ns-1)];
    Out1 = fft(y1a, Ns);
    a1 = abs(Out1); 
    a1(Ns/2 : end) = 0;
    %figure; plot(a1)
    [max1_abs(i), m_pos1(i)] = max(a1);
    p1 = phase(Out1(m_pos1(i)));
    f1 = m_pos1(i)/Ns*fs;
%     if i == 181
%         figure; plot(a1)
%     end
s2 = x_m_c(2*Ns+1:4*Ns-1);
    y2 = s2 .* win2;
    y2a = y2(Ns:end) + [0; y2(1:Ns-1)];
    Out2 = fft(y2a, Ns);
    a2 = abs(Out2);
    a2(Ns/2 : end) = 0;
    %figure; plot(a2)
    [max2_abs(i), m_pos2(i)] = max(a2);    
     p2 = phase(Out2(m_pos2(i)));
     f2 = m_pos2(i)/Ns*fs;
%-------俩个序列都做apfft

%-------第一个序列做fft
% s1 = x_m_c(Ns : 2*Ns-1);
%     y1 = s1 .* win_1;
%     Out1 = fft(y1, Ns);
%     a1 = abs(Out1); 
%     a1(Ns/2 : end) = 0;
%     %figure; plot(a1)
%     [max1_abs(i), m_pos1(i)] = max(a1);
%     p1 = phase(Out1(m_pos1(i)));
%     f1 = m_pos1(i)/Ns*fs;
%-------第一个序列做fft
%-------第二个序列做apfft
% s2 = x_m_c(1:2*Ns-1);
%     y2 = s2 .* win_2;
%     y2a = y2(Ns:end) + [0; y2(1:Ns-1)];
%     Out2 = fft(y2a, Ns);
%     a2 = abs(Out2);
%     a2(Ns/2 : end) = 0;
%     %figure; plot(a2)
%     [max2_abs(i), m_pos2(i)] = max(a2);    
%      p2 = phase(Out2(m_pos2(i)));
%-------第二个序列做apfft
%-------双全相位法的频率校正
    %fhi(i) = f1-f2;
    phi(i) = p2-p1;
    if phi(i) > pi
        phi(i) = phi(i) - 2*pi;
    else if phi(i) < -pi
            phi(i) = phi(i) + 2*pi;
        end
    end
    fi(i) = round(2*Ns/fs*f1)*2*pi;
    hi(i) = fi(i) + phi(i);
    f_fine1(i) = hi(i)/(2*Ns)/2/pi*fs;
%-------双全相位法的频率校正

%-------fft-apfft的频率校正
%     phi(i) = p1-p2;
%     if phi(i) > pi
%         phi(i) = phi(i) - 2*pi;
%     else if phi(i) < -pi
%             phi(i) = phi(i) + 2*pi;
%         end
%     end
%     hi(i) = phi(i) / (Ns-1) * 2;   %频偏--角频率
%     fi(i) = f1;  %--加窗传统fft观测到的频率
%     f_fine1(i) = hi(i)/pi/2*fs + f1;
%-------fft-apfft的频率校正

end

mean(f_fine1)
std(f_fine1)
maxo = max(f_fine1) - min(f_fine1)
figure; plot(m_pos1)
figure; plot(f_fine1)
figure; plot(hi)
figure; plot(phi)
figure; plot(fi)
