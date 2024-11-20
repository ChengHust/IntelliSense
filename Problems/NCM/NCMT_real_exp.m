classdef NCMT_real_exp < PROBLEM

    properties
        Data; % Measured data
        O; % Optimal decision vector
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
             [data_m,data_ang]= obj.ParameterSet([],[]);
             obj.lower = [9,9,9,0,0,0,0,0,0,0,2*pi/3,4*pi/3];
             obj.upper = [12,12,12,2*pi,2*pi,2*pi,25,25,25,5*pi/180,125*pi/180,245*pi/180];
             obj.D = 12;
             obj.M = 1;
             obj.N = 100;
             obj.Data = [data_m;data_ang];
             obj.encoding = ones(1,obj.D);
        end
        
        function PopObj= CalObj(obj,PopDec)
            n      = 3;           
            rs     = 39.5e-3;
            N_cal  = 1001;
            u0     = 4*pi*10e-7;
            [N,~]  = size(PopDec);
            LocDec = PopDec(:,1:6);
            CurDec = PopDec(:,7:end);
            LocSen_p  = rs*ones(N, 3*N_cal);
            LocSen_salph = repmat(linspace(pi/2, 5*pi/2, N_cal) , N , 3);
            LocDec_p    = repelem(LocDec(:,1:3),1,N_cal)/1000;
            LocDec_aplh = repelem(LocDec(:,4:6),1,N_cal);
            I=repelem([CurDec(:,1).*exp(1i*CurDec(:,1+n)),CurDec(:,2).*exp(1i*CurDec(:,2+n)),CurDec(:,3).*exp(1i*CurDec(:,3+n))],1,N_cal);
            B_t=u0*I/(2*pi).*((LocSen_p-LocDec_p.*cos(LocSen_salph-LocDec_aplh))./(LocSen_p.^2+LocDec_p.^2-2*LocSen_p.*LocDec_p.*cos(LocSen_salph-LocDec_aplh)));
            B_x=(B_t(:,1:N_cal)+B_t(:,N_cal+1:2*N_cal)+B_t(:,2*N_cal+1:3*N_cal))/10;
            PopObj = cal_fit(B_x, obj.Data(1,:),obj.Data(2,:));
        end
    end 
end

function fitness = cal_fit(data_c, data_m,data_ang)
[data_c_row_N,data_c_column_N] = size(data_c);
[~,data_m_column_N] = size(data_m);
[~,data_ang_column_N] = size(data_ang);
data_m_r = signal_rec(data_m,data_c_column_N,data_m_column_N);
n=17;
n_a = floor(linspace(1,1001, n + 1));
data_c_sam_ang =angle(data_c(:,n_a(1:n)));

e_fz = abs(repmat(abs(data_m_r),data_c_row_N,1)- abs(data_c)) ./ abs(data_c);
e_xw = abs(repmat(abs(data_ang),data_c_row_N,1)- abs(data_c_sam_ang)) ./ abs(data_c_sam_ang);

fitness_fz= sum(e_fz,2)/ data_c_column_N;
fitness_xw= sum(e_xw,2)/data_ang_column_N;
fitness = fitness_fz+fitness_xw*10e-5;
end