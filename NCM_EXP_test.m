function NCM_EXP_test()
clear,clc;format compact
load real_exp_data.mat
alg = {@GA,@PSO,@DE,@MVPA,@ECPO,@IMODE,@IMODE_Fast_real};

for i=1:5
    for j=1:7
        platemo('algorithm',alg {j},'N',100,'M',1,'problem',{@NCMT_real_exp,data_abs(i,:),data_angle(i,:)},'maxFE',40000,'save',20);
    end
end
end


