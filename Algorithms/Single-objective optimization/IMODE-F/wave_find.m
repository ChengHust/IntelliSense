function Ang= wave_find(data) %��������data������ֱ�Ϊ����Ƕȡ�����λ�á�����λ��
%�ź���ɢ��
data_exp=[data,data];
IndMax=find(diff(sign(diff(abs(data_exp))))<0)+1;
IndMax=IndMax(1:3);
Ang=IndMax.*2*pi/length(data); 
end