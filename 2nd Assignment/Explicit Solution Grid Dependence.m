%%%%%%%%%%%%% CFD 2ND ASSIGMENT EXPLICIT SCHEME GRID DEPENDENCE VICHAS PAVLOS AEM: 6663 %%%%%%%%%%%%%%%
clear
clc
close all

Reynolds = [1,100,500];
YXpoints = [30 ,50,80,100,120 ];

for k = 1:5
    Xpoints = YXpoints(k);
for g = 1:3

    Re = Reynolds(g);

U = 1;        %%% ΤΑΧΥΤΗΤΑ ΑΝΩ ΠΛΑΚΑΣ
V = 1; %%% KINEMATIC VISCOCITY

Lx = 1;  %%% ΜΗΚΟΣ ΚΟΙΛΩΜΑΤΟΣ m
Ly = 1;  %%% ΥΨΟΣ ΚΟΙΛΩΜΑΤΟΣ  m 

%Xpoints = 160; %%%ΣΗΜΕΙΑ ΣΤΟΝ ΑΞΟΝΑ Χ (j points) ΣΤΗΛΕΣ
%Ypoints = 160; %%%ΣΗΜΕΙΑ ΣΤΟΝ ΑΞΟΝΑ Y (i points) ΓΡΑΜΜΕΣ

%Re = (Lx*U)/V;          %%% REYNOLDS NUMBER
Ypoints = Xpoints;

if Re==1
dt    = 0.00001;
Relax = 1.2;        %%% Xpoints=Ypoints=100 & 5000 iterations
iter = 5000;


elseif Re==100
dt    = 0.001;
Relax = 0.9;        %%% Xpoints=Ypoints=100 & 15000 iterations
iter = 15000;


elseif Re == 500
dt    = 0.001;    
Relax = 1.5;          %%% Xpoints=Ypoints=160 & 100000 iterations
iter = 180000;

end


StepX = (Lx/(Xpoints-1));   %%% ΒΗΜΑ ΣΤΗΝ ΚΑΤΕΥΘΘΥΝΣΗ Χ (j) ΣΤΗΛΕΣ
StepY = (Ly/(Ypoints-1));   %%% ΒΗΜΑ ΣΤΗΝ ΚΑΤΕΥΘΘΥΝΣΗ Y (i) ΓΡΑΜΜΕΣ

%%%% INITIAL GRIDS %%%%
Vort = zeros(Ypoints,Xpoints);
Psi  = zeros(Ypoints,Xpoints);
u    = zeros(Ypoints,Xpoints);
v    = zeros(Ypoints,Xpoints);


%%%%% BOUNDARY CONDITIONS %%%%%
% FOR STREAM FUNCTION
Psi(:,1) = 0;
Psi(:,end) = 0;
Psi(1,:) = 0;
Psi(end,:) = 0;

% FOR VORTICITY
Vort(:,1) = 2*(Psi(:,1)-Psi(:,2))/(StepX^2);                    %%% LEFT WALL
Vort(:,end) = 2*(Psi(:,end)-Psi(:,end-1))/(StepX^2);            %%% RIGHT WALL
Vort(end,:) = 2*(Psi(end,:)-Psi(end-1,:))/(StepY^2);            %%% BOTTOM WALL
Vort(1,:) = (2*(Psi(1,:)-Psi(2,:))/(StepY^2)) - (2*U/StepY);    %%% TOP WALL

%%% FOR RENEW MATRIX %%%
Vort_new = Vort;
Psi_new  = Psi;
maxiterations = dt*iter;     %%% max iterations
timer = 1;                   %%% just a timer
Courant = (dt*U)/StepX; %%% STABILLITY ANALYSIS


tic
%%% CALCULATIONS %%%

for time = 0:dt:maxiterations
     
    Vort_new = Vort + Relax*(Vort_new - Vort); %%%% RELAXATION
    
    Vort = Vort_new;
    Psi  = Psi_new;
    
    %%% RENEW BOUNDARY CONDITIONS %%%
    Vort(:,1) = 2*(Psi(:,1)-Psi(:,2))/(StepX^2);                    %%% LEFT WALL
    Vort(:,end) = 2*(Psi(:,end)-Psi(:,end-1))/(StepX^2);            %%% RIGHT WALL
    Vort(end,:) = 2*(Psi(end,:)-Psi(end-1,:))/(StepY^2);            %%% BOTTOM WALL
    Vort(1,:) = (2*(Psi(1,:)-Psi(2,:))/(StepY^2)) - (2*U/StepY);    %%% TOP WALL
    
    
    for j = 2:Xpoints-1 %%% ΣΤΗΛΕΣ
        
        for i = 2:Ypoints-1  %%% ΓΡΑΜΜΕΣ
         
            %%% VORTICITY EQUATION %%%
            d2Vortdx2 = (Vort(i,j+1) - 2*Vort(i,j) + Vort(i,j-1))/(StepX^2);
            d2Vortdy2 = (Vort(i+1,j) - 2*Vort(i,j) + Vort(i-1,j))/(StepY^2);
            dPsidy_dVortdx = 0.25*((Psi(i+1,j)-Psi(i-1,j))/StepY)*((Vort(i,j+1)-Vort(i,j-1))/(StepX));
            dPsidx_dVortdy = 0.25*((Psi(i,j+1)-Psi(i,j-1))/StepX)*((Vort(i+1,j)-Vort(i-1,j))/(StepY));
        
            Vort_new(i,j) = Vort(i,j) + dt*((d2Vortdx2+d2Vortdy2)*(1/Re) - dPsidy_dVortdx +   dPsidx_dVortdy);  
        
            %%% PSI EQUATION %%%
            Psi_new(i,j) = (1/4)*(Vort_new(i,j)*(StepX^2) + (Psi(i,j+1)+Psi(i,j-1)) +  (Psi(i+1,j)+Psi(i-1,j)));
            
            
        end    
    
    end
   
%     if time==dt
%         Vort_initial = max(Vort_new);
%     end    
%     
    %%% Check for Converge %%%
    if time>dt*10
    ERROR(timer) = max(max(Vort_new - Vort));
    if ERROR(timer) < 1e-4
        fprintf('man μου οι επαναληψεις είναι %d με βήμα %d\n', time/dt, dt);
        break
    end
    timer=timer+1;
    end
    
    su(timer) = sum(sum(Vort_new));
    
    
    
end
toc


Array_Vort(k,g) = Vort_new(Xpoints/2,Xpoints/2);  %%% Store Vorticity values
Array_Psi(k,g) = Psi_new(Xpoints/2,Xpoints/2);  %%% Store Vorticity values

end
end

%%%% PLOT GRID DEPENDENCE ANALYSIS %%%%
%%% Vorticity values %%%
figure
plot([30,50,80,100,120],Array_Vort(:,1),'r','LineWidth',2)
hold on
plot([30,50,80,100,120],Array_Vort(:,2),'b','LineWidth',2)
hold on
plot([30,50,80,100,120],Array_Vort(:,3),'k','LineWidth',2)
xlabel('Xpoints of Square Grid (Xpoints = Ypoints)')
ylabel('Vorticity values')
title('Explicit scheme Grid dependence (Στροβιλότητα)')
legend('Re = 1','Re = 100','Re = 500')
grid on

%%% Vorticity values %%%
figure
plot([30,50,80,100,120],Array_Psi(:,1),'r','LineWidth',2)
hold on
plot([30,50,80,100,120],Array_Psi(:,2),'b','LineWidth',2)
hold on
plot([30,50,80,100,120],Array_Psi(:,3),'k','LineWidth',2)
xlabel('Xpoints of Square Grid (Xpoints = Ypoints)')
ylabel('Psi function values')
title('Explicit scheme Grid dependence (Ροική Συνάρτηση)')
legend('Re = 1','Re = 100','Re = 500')
grid on
