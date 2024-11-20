function Ang= wave_find(data) %输入数据data；输出分别为导体角度、波峰位置、波谷位置
%信号离散化
data_exp=[data,data];
IndMax=find(diff(sign(diff(abs(data_exp))))<0)+1;
IndMax=IndMax(1:3);
Ang=IndMax.*2*pi/length(data); 
end