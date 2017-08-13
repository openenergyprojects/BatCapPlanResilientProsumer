clear all;
clc;

%input data
% initial_capacity_battery=0; % KWh
% minimum_capacity_battery=0; %KWh
% maximum_capacity_battery=8; %KWh
% 
% max_batt_discharge=3;%KW
% max_batt_charge=3; %KW
% 
% battery_effic_disch=0.95; %it cant be 1. set 0.99
% battery_effic_charge=0.95; %it cant be 1. set 0.99

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TESLA BATTERY
%input data
initial_capacity_battery=0; % KWh
minimum_capacity_battery=0; %KWh
maximum_capacity_battery=13.5; %KWh

max_batt_discharge=5;%KW
max_batt_charge=5; %KW

battery_effic_disch=0.95; %it cant be 1. set 0.99
battery_effic_charge=0.95; %it cant be 1. set 0.99

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%iNVERTER

inverter_size=5;%KW

%discretization of the converter efficiency
efficiencies=[0.93 0.95 0.96 0.975];
bounds_efficiency_temp=[0 10 20 30]; %percentage




%Solver option: 1:Gurobi, 2:Matlab
solver_option=1;

% path_data='micro_grid_data.xlsx';
a=importdata('micro_grid_dataIC.csv'); %much faster (open at once the file) than reading with xlsread (open file each time you read a new variable :D)
%end of input data
%% ***************************************************************************** 
%data processing
bet=bounds_efficiency_temp;
bounds_efficiency=[bet(1)*inverter_size/100 bet(2)*inverter_size/100 bet(3)*inverter_size/100 bet(4)*inverter_size/100]; %KW
%bounds_efficiency=[0 4 5 6];

% read the input data
% tic
% Load=xlsread(path_data,'A1:A8'); %Load (KW) per each quarter hour
% Load=Load';
% PV=xlsread(path_data,'B1:B8'); %PV generation (KW) per each quarter hour
% PV=PV';
% c=xlsread(path_data,'C1:C8'); %Price of buy 1KWh per each quarter hour
% c=c';
% k=xlsread(path_data,'D1:D8'); %Price of sell 1KWh per each quarter hour
% k=k';
% c=c/4; %Price of buy 1KW continually for 15 mins
% k=k/4; %Price of sell 1KW continually for 15 mins
%  
% 
% read_time=toc



% tic
% Load=xlsread(path_data,'N1:N96'); %Load (KW) per each quarter hour
% Load=Load';
% PV=xlsread(path_data,'O1:O96'); %PV generation (KW) per each quarter hour
% PV=PV';
% c=xlsread(path_data,'P1:P96'); %Price of buy 1KWh per each quarter hour
% c=c';
% k=xlsread(path_data,'Q1:Q96'); %Price of sell 1KWh per each quarter hour
% k=k';
% c=c/4; %Price of buy 1KW continually for 15 mins
% k=k/4; %Price of sell 1KW continually for 15 mins
%  
% 
% read_time=toc

%%%%%%%%%%%%%%%%%%%%
tic
Load=xlsread(path_data,'U1:U192'); %Load (KW) per each quarter hour
Load=Load';
PV=xlsread(path_data,'V1:V192'); %PV generation (KW) per each quarter hour
PV=PV';
c=xlsread(path_data,'W1:W192'); %Price of buy 1KWh per each quarter hour
c=c';
k=xlsread(path_data,'X1:X192'); %Price of sell 1KWh per each quarter hour
k=k';
c=c/4; %Price of buy 1KW continually for 15 mins
k=k/4; %Price of sell 1KW continually for 15 mins
 
read_time=toc





num_descrit=length(efficiencies);

%end of input data
%% *****************************************************************************    
%prepare the matrixes
icb=initial_capacity_battery;
mincb=minimum_capacity_battery;
maxcb=maximum_capacity_battery;


num_of_hours=length(c);

%1D zero element matrix (time intervals)
zeros_1D=zeros(1,num_of_hours);

%2D zero element matrix (time intervals*time intervals)
zeros_2D=zeros(num_of_hours,num_of_hours);

%Create the positive diagonal matrix with 1.
v = ones(1,num_of_hours);
Diag1_pos = diag(v);

%Create the negative diagonal matrix with -1.
v = ones(1,num_of_hours);
v=v.*-1;
Diag1_neg =diag(v);

%Create the negative diagonal battery efficiency.
Diag_neg_disch_eff=Diag1_neg*battery_effic_disch;  %is not implemented yet
Diag_pos_charge_eff=Diag1_pos*(1/battery_effic_charge);
%Create the  diagonal matrix OF THE first efficiency
v = ones(1,num_of_hours);
v=v.*efficiencies(1);
Diag_eff1_pos =diag(v);
Diag_eff1_neg=-Diag_eff1_pos;

%Create the  diagonal matrix OF THE second efficiency
v = ones(1,num_of_hours);
v=v.*efficiencies(2);
Diag_eff2_pos =diag(v);
Diag_eff2_neg=-Diag_eff2_pos;

%Create the  diagonal matrix OF THE third efficiency
v = ones(1,num_of_hours);
v=v.*efficiencies(3);
Diag_eff3_pos =diag(v);
Diag_eff3_neg=-Diag_eff3_pos;

%Create the  diagonal matrix OF THE fourth efficiency
v = ones(1,num_of_hours);
v=v.*efficiencies(4);
Diag_eff4_pos =diag(v);
Diag_eff4_neg=-Diag_eff4_pos;

%create lower triangular matrix positive and negative
v=tril(ones(num_of_hours,num_of_hours),-1);
pos_triangular=v+Diag1_pos;
neg_triangular=-pos_triangular;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos_triangular=pos_triangular/4; %for quarter hour
neg_triangular=neg_triangular/4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create the  diagonal matrix OF THE bounds efficiency 
v = ones(1,num_of_hours);
v=v.*bounds_efficiency(1);
Diag_bounds_eff1 =diag(v);

%Create the  diagonal matrix OF THE bounds efficiency 
v = ones(1,num_of_hours);
v=v.*bounds_efficiency(2);
Diag_bounds_eff2 =diag(v);

%Create the  diagonal matrix OF THE bounds efficiency 
v = ones(1,num_of_hours);
v=v.*bounds_efficiency(3);
Diag_bounds_eff3 =diag(v);

%Create the  diagonal matrix OF THE bounds efficiency 
v = ones(1,num_of_hours);
v=v.*bounds_efficiency(4);
Diag_bounds_eff4 =diag(v);

%create the upper bound value of y1 and y2
upper_y1=max_batt_discharge+max(PV);
upper_y2=max_batt_charge*(1/battery_effic_charge);

%Create the  diagonal matrix of upper bound value of y1 
v = ones(1,num_of_hours);
v=v.*upper_y1; %very big value for ours values
Diag_upper_pos_y1 =diag(v);
Diag_upper_neg_y1 =-Diag_upper_pos_y1;

%Create the  diagonal matrix of upper bound value of y2 
v = ones(1,num_of_hours);
v=v.*upper_y2; %very big value for ours values
Diag_upper_pos_y2 =diag(v);
Diag_upper_neg_y2 =-Diag_upper_pos_y2;
%% *****************************************************************************    
 %Objective Function
%
 % F=[pgrid_pos pgrid_neg bat PDC y1:y4 B1:B4];
 F=[c k*(-1) zeros(1,num_of_hours*4) zeros(1,num_of_hours*4*num_descrit)];

intcon = 6*num_of_hours+num_of_hours*2*num_descrit+1:length(F); %binary variables
% %BOUNDS


lb=[zeros(1,2*num_of_hours) zeros(1,2*num_of_hours)  zeros(1,2*num_of_hours)  zeros(1,num_of_hours*4*num_descrit) ];
ub=[ones(1,2*num_of_hours)*(inf)  ones(1,num_of_hours)*(max_batt_discharge) ones(1,num_of_hours)*(max_batt_charge) ones(1,2*num_of_hours)*(inf) ones(1,num_of_hours*4*num_descrit)*(inf) ];

 
 %Equalities Constraints
Aeq=[Diag1_pos,Diag1_neg,zeros_2D,zeros_2D,zeros_2D,Diag1_neg,Diag1_pos,Diag1_pos,Diag1_pos,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D;
    zeros_2D,zeros_2D,Diag_neg_disch_eff,Diag_pos_charge_eff,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_neg,Diag1_neg,Diag1_neg,Diag1_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D ;
    zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,Diag1_pos,Diag1_pos,Diag1_pos,Diag1_pos,Diag1_pos,Diag1_pos,Diag1_pos;
    ];

 beq=[Load,PV,ones(1,num_of_hours)];

 %Inequalities Constraints

 A=[zeros_2D,zeros_2D,pos_triangular,neg_triangular,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D;%5constraint
     zeros_2D,zeros_2D,neg_triangular,pos_triangular,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D ;%5constraint
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_bounds_eff2,Diag_bounds_eff3,Diag_bounds_eff4,zeros_2D,zeros_2D,zeros_2D,zeros_2D;%7constraint
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_bounds_eff2,Diag_bounds_eff3,Diag_bounds_eff4;%8constraint
    %
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_neg_y1,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %9constraint Y1.1,B1
    zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_neg_y1,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %9constraint Y1.2,B2
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_neg_y1,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %9constraint Y1.3,B3
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_neg_y1,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %9constraint Y1.4,B4
     %
      zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_neg_y2,zeros_2D,zeros_2D,zeros_2D; %9constraint Y2.1,Z1
    zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_neg_y2,zeros_2D,zeros_2D; %9constraint Y2.2,Z2
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_neg_y2,zeros_2D; %9constraint Y2.3,Z3
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_neg_y2; %9constraint Y2.4,Z4
     %
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff1_neg,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %10constraint Y1.1,X1
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff2_neg,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %10constraint Y1.2,X1
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff3_neg,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %10constraint Y1.3,X1
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff4_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %10constraint Y1.4,X1
     %
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff1_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %10constraint Y2.1,X2
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff2_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %10constraint Y2.2,X2
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff3_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %10constraint Y2.3,X2
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff4_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %10constraint Y2.4,X2
     %
      %
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff1_pos,zeros_2D,Diag1_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_pos_y1,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %11constraint Y1.1,X1,B1
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff2_pos,zeros_2D,zeros_2D,Diag1_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_pos_y1,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %11constraint Y1.2,X1,B1
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff3_pos,zeros_2D,zeros_2D,zeros_2D,Diag1_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_pos_y1,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %11constraint Y1.3,X1,B1
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff4_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_pos_y1,zeros_2D,zeros_2D,zeros_2D,zeros_2D; %11constraint Y1.4,X1,B1
     %
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_pos_y2,zeros_2D,zeros_2D,zeros_2D; %11constraint Y2.1,X2,B2
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff2_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_pos_y2,zeros_2D,zeros_2D; %11constraint Y2.2,X2,B2
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff3_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_pos_y2,zeros_2D; %11constraint Y2.3,X2,B2
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_eff4_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_neg,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_upper_pos_y2; %11constraint Y2.4,X2,B2
    
      zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_bounds_eff2*(-1),Diag_bounds_eff3*(-1),Diag_bounds_eff4*(-1),Diag1_neg*inverter_size,zeros_2D,zeros_2D,zeros_2D,zeros_2D;%12constraint
     zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag1_pos,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,zeros_2D,Diag_bounds_eff2*(-1),Diag_bounds_eff3*(-1),Diag_bounds_eff4*(-1),Diag1_neg*inverter_size;%13constraint
     ];
     
 
 b=[ones(1,num_of_hours)*(icb-mincb),ones(1,num_of_hours)*(maxcb-icb),zeros(1,2*num_of_hours),zeros(1,2*8*num_of_hours),ones(1,4*num_of_hours)*upper_y1,ones(1,4*num_of_hours)*upper_y2,zeros(1,2*num_of_hours)];
 
% %% *****************************************************************************    
%% *****************************************************************************    
%Call the solver
if solver_option==2
%Call Matlab MILP
    tic;
    [x,fval] = intlinprog(F,intcon,A,b,Aeq,beq,lb,ub);
    toc;

elseif   solver_option==1
%Call Gurobi MILP
    path_Gurobi='gurobi';
     addpath(path_Gurobi) %if you choose gurobi for dispatch, then run gurobi_setup
            savepath();
    f=F;

    nargin=8;
    if nargin < 4
        error('intlinprog(f, intcon, A, b)')
    end

    if nargin > 8
        error('intlinprog(f, intcon, A, b, Aeq, beq, lb, ub)');
    end

    if ~isempty(A)
        n = size(A, 2);
    elseif nargin > 5 && ~isempty(Aeq)
        n = size(Aeq, 2);
    else
        error('No linear constraints specified')
    end

    if ~issparse(A)
        AA = sparse(A);
    end

    if nargin > 4 && ~issparse(Aeq)
        Aeqq = sparse(Aeq);
    end

    model.obj = f;
    model.vtype = repmat('C', n, 1);
    model.vtype(intcon) = 'I';

    if nargin < 5
        model.A = AA;
        model.rhs = b;
        model.sense = '<';
    else
        model.A = [AA; Aeqq];
        model.rhs = [b'; beq'];
        model.sense = [repmat('<', size(A,1), 1); repmat('=', size(Aeq,1), 1)];
    end

    if nargin < 7
        model.lb = -inf(n,1);
    else
        model.lb = lb;
    end

    if nargin == 8
       model.ub = ub;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  params.MIPGap=1e-3;
   %  params.FeasibilityTol=1e-4;
   % params.IntFeasTol=1e-3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    params.outputflag = 1;
    tic
    result = gurobi(model, params);
    toc

    if strcmp(result.status, 'OPTIMAL')
        exitflag = 1;
    elseif strcmp(result.status, 'INTERRUPTED')
        if isfield(result, 'x')
            exitflag = 2;
        else
            exitflag = 0;
        end
    elseif strcmp(result.status, 'INF_OR_UNBD')
        params.dualreductions = 0;
        result = gurobi(model, params);
        if strcmp(result.status, 'INFEASIBLE')
            exitflag = -2;
        elseif strcmp(result.status, 'UNBOUNDED')
            exitflag = -3;
        else
            exitflag = nan;
        end
    else
        exitflag = nan;
    end


    if isfield(result, 'x')
        x = result.x;
    else
        x = nan(n,1);
    end

    if isfield(result, 'objval')
        fval = result.objval;
    else
        fval = nan;
    end

end; %if

%% *****************************************************************************    
%Total COST without the battery
% Buy energy= Load -PV (if it is negative then you sell energy)

for i=1:num_of_hours
    if PV(i)< bounds_efficiency(2)
    PV_NET(i)=PV(i)*efficiencies(1);
    elseif PV(i)>= bounds_efficiency(2) && PV(i)< bounds_efficiency(3)
    PV_NET(i)=PV(i)*efficiencies(2);
    elseif PV(i)>= bounds_efficiency(3) && PV(i)< bounds_efficiency(4)
    PV_NET(i)=PV(i)*efficiencies(3);
    else    
    PV_NET(i)=PV(i)*efficiencies(4);
    end;
end;
Cost_without_battery=0;
for i=1:num_of_hours
    NET_DEMAND(i)=Load(i)-PV_NET(i);
    if NET_DEMAND(i)>=0
       Cost_without_battery=Cost_without_battery+NET_DEMAND(i)*c(i);
    else    
       Cost_without_battery=Cost_without_battery+NET_DEMAND(i)*k(i);
    end;    
end;
fprintf('Cost_without_battery=%d',Cost_without_battery)
Cost_without_battery
%% *****************************************************************************    

Pbuy_interval=x(1:num_of_hours);
Psell_interval=x(num_of_hours+1:2*num_of_hours);
Pbattery_inter_disch=x(2*num_of_hours+1:3*num_of_hours);
Pbattery_inter_charge=x(3*num_of_hours+1:4*num_of_hours);
Px1_interval=x(4*num_of_hours+1:5*num_of_hours);
Px2_interval=x(5*num_of_hours+1:6*num_of_hours);
Py11_interval=x(6*num_of_hours+1:7*num_of_hours);
Py12_interval=x(7*num_of_hours+1:8*num_of_hours);
Py13_interval=x(8*num_of_hours+1:9*num_of_hours);
Py14_interval=x(9*num_of_hours+1:10*num_of_hours);
Py21_interval=x(10*num_of_hours+1:11*num_of_hours);
Py22_interval=x(11*num_of_hours+1:12*num_of_hours);
Py23_interval=x(12*num_of_hours+1:13*num_of_hours);
Py24_interval=x(13*num_of_hours+1:14*num_of_hours);
B1_interval=x(14*num_of_hours+1:15*num_of_hours);
B2_interval=x(15*num_of_hours+1:16*num_of_hours);
B3_interval=x(16*num_of_hours+1:17*num_of_hours);
B4_interval=x(17*num_of_hours+1:18*num_of_hours);
Z1_interval=x(18*num_of_hours+1:19*num_of_hours);
Z2_interval=x(19*num_of_hours+1:20*num_of_hours);
Z3_interval=x(20*num_of_hours+1:21*num_of_hours);
Z4_interval=x(21*num_of_hours+1:22*num_of_hours);



%TOTAL COST OF OPERATION
Total_Cost_of_Operation=c*Pbuy_interval-k*Psell_interval %cents


 
figure(1);
subplot(2,2,3)
stem(Pbuy_interval-Psell_interval,'BaseValue',0)
title('Grid Buy/Sell Energy')
xlabel('Time intervals')
ylabel(' Buy/Sell(KW)')

subplot(2,2,1)
stem(Load)
title(' Load')
xlabel('Time intervals')
ylabel('(KW)')

subplot(2,2,2)
stem(c,'BaseValue',0)
hold on;
stem(-k,'BaseValue',0)
title(' Buy/Sell Cost')
xlabel('Time intervals')
ylabel(' Buy/Sell(cents/KW for 15min)')

subplot(2,2,4)
stem(Pbattery_inter_disch-Pbattery_inter_charge,'BaseValue',0)
title('Charge/Discharge Energy')
xlabel('Time intervals')
ylabel(' Charge/Discharge(KW)')
 

%TO DIA 4 EN GIA QUARTER THS WRAS
battery_power_per_hour(1)=icb-Pbattery_inter_disch(1)/4+Pbattery_inter_charge(1)/4;
for i=2:num_of_hours
battery_power_per_hour(i)=battery_power_per_hour(i-1)-Pbattery_inter_disch(i)/4+Pbattery_inter_charge(i)/4;
end;

figure(2);
subplot(2,1,1)
stem(battery_power_per_hour)
title('Capacity level of the battery for each quarter of hour')
xlabel('Time intervals')
ylabel('(KWh)')

subplot(2,1,2)
stem(PV)
title('PV Generation for each quarter of hour')
xlabel('Time intervals')
ylabel('(KW)')

x_time=[0:1/4:48];
figure(4)
plot(x_time(1:end-1),Pbuy_interval, '.-b',...
     x_time(1:end-1),Load, '.-r',...
     x_time(1:end-1),PV, '.-y',...
     x_time(1:end-1),Pbattery_inter_disch, '.-k', ...
     x_time(1:end-1),-Pbattery_inter_charge, '.-g');
title('Scheduling of Grid Energy')
xlabel('Time (h)')
ylabel ('Power (kW)')
legend('Pgrid', 'Pload', 'Ppv', 'Pbat-Disch', 'Pbat-Charge');