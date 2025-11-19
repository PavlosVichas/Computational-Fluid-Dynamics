%%%%%%%%%%%%% CFD 2ND ASSIGMENT EXPLICIT SCHEME VICHAS PAVLOS AEM: 6663 %%%%%%%%%%%%%%%
clear
clc

U = 1;        %%% ΤΑΧΥΤΗΤΑ ΑΝΩ ΠΛΑΚΑΣ %%% - (2*U/StepY) for left + (2*U/StepY) for right
V = 0.002;    %%% KINEMATIC VISCOCITY


Lx = 1;  %%% ΜΗΚΟΣ ΚΟΙΛΩΜΑΤΟΣ m
Ly = 1;  %%% ΥΨΟΣ ΚΟΙΛΩΜΑΤΟΣ  m 

%Xpoints =50; %%%ΣΗΜΕΙΑ ΣΤΟΝ ΑΞΟΝΑ Χ (j points) ΣΤΗΛΕΣ
%Ypoints = 50; %%%ΣΗΜΕΙΑ ΣΤΟΝ ΑΞΟΝΑ Y (i points) ΓΡΑΜΜΕΣ

%Re = (Lx*U)/V;          %%% REYNOLDS NUMBER
Reynolds = [1,100];
YXpoints=[30,40,50,60,100];

for k = 1:4
    
Xpoints = YXpoints(k);    
  Ypoints = Xpoints;      
for g = 1:2
Re = Reynolds(g);
   
if Re==1
dt = 0.001;
%Xpoints = 50;
%Ypoints = 50;
max_iterations = 8000;       
Relax = 0.8;
end
if Re == 100
dt = 0.001;                   
max_iterations = 800;       
Relax = 0.9;
% Xpoints = 50;
% Ypoints = 50;    
end

StepX = (Lx/(Xpoints-1));   %%% ΒΗΜΑ ΣΤΗΝ ΚΑΤΕΥΘΘΥΝΣΗ Χ (j) ΣΤΗΛΕΣ
StepY = (Ly/(Ypoints-1));   %%% ΒΗΜΑ ΣΤΗΝ ΚΑΤΕΥΘΘΥΝΣΗ Y (i) ΓΡΑΜΜΕΣ

%%%% INITIAL GRIDS %%%%
Vort = zeros(Ypoints,Xpoints);
Psi  = zeros(Ypoints,Xpoints);
u    = zeros(Ypoints,Xpoints);
v    = zeros(Ypoints,Xpoints);
K_half = zeros(Ypoints-2,1);
Aa = zeros(Ypoints-2,1);
Bb = zeros(Ypoints-3,1);
K_full = zeros(1,Xpoints-2);
K_half_Psi = zeros(Ypoints-2,1);
VORT_HALF = zeros(Ypoints,Xpoints-2); %%% -2 Γιατί μετά κουμπώνω τις οριακές συνθήκες
VORT_FULL = zeros(Ypoints-2,Xpoints);

%%%%% BOUNDARY CONDITIONS %%%%%
% FOR STREAM FUNCTION
Psi(:,1) = 0;
Psi(:,end) = 0;
Psi(1,:) = 0;
Psi(end,:) = 0;
Psi_new=Psi;

timer = 1;                   %%% just a timer
tic

%%% IMPLICIT EULER SCHEME %%%

for time = 0:dt:max_iterations
    
    Vort_previous = Vort;  %%% Tool for converge
    Psi_previous   = Psi;
    %Psi=Psi_new;   
    % FOR VORTICITY

Vort(:,1) = 2*(Psi(:,1)-Psi(:,2))/(StepX^2);                    %%% LEFT WALL
Vort(:,end) = 2*(Psi(:,end)-Psi(:,end-1))/(StepX^2);            %%% RIGHT WALL
Vort(end,:) = 2*(Psi(end,:)-Psi(end-1,:))/(StepY^2);            %%% BOTTOM WALL
Vort(1,:) = (2*(Psi(1,:)-Psi(2,:))/(StepY^2)) - (2*U/StepY);    %%% TOP WALL 

    
    
    %%%%%% HALF TIME n +1/2 MARCHING IN Y-DIRECTION  %%%%%%
    
    %%% Vort matrix coefficients constants %%%
    C = (dt/2)*(1/Re)*(1/(StepY^2));
    D = (dt/2)*(1/Re)*(1/(StepX^2));
    
    
    %%% K matrix calculations for Vort %%%
for j = 2:Xpoints-1
    
for i = 2:Ypoints -1   
       
       %%% Variable Coeficients %%%
       A(i-1) = (dt/2)*(1/(StepX*StepY))*(Psi(i+1,j) - Psi(i,j)); 
       B(i-1) = (dt/2)*(1/(StepX*StepY))*(Psi(i,j+1) - Psi(i,j)); 
          
    
      if i == 2               %%% Για πρωτη γραμμη ξέρω vorticity
      K_half(i-1,1) = D*Vort(i,j-1) + (1+A(i-1)-2*D)*Vort(i,j) + (-A(i-1)+D)*Vort(i,j+1) + C*Vort(i-1,j);
      elseif i == Ypoints-1   %%% Για τελευταία γραμμη ξέρω vorticity
      K_half(i-1,1) = D*Vort(i,j-1) + (1+A(i-1)-2*D)*Vort(i,j) + (-A(i-1)+D)*Vort(i,j+1) + (C+B(i-1))*Vort(i+1,j);
      else
      K_half(i-1,1) = D*Vort(i,j-1) + (1+A(i-1)-2*D)*Vort(i,j) + (-A(i-1)+D)*Vort(i,j+1);    
      end
end
      

          
      
      %%% Construct diagonals %%%
       Coef1 = -C*ones(Ypoints - 3 , 1);
       Coef2 = (1+B+2*C)'.*ones(Ypoints - 2 , 1);
       Coef3 = (-B(1:end-1)-C)'.*ones(Ypoints - 3 , 1);
      
       AA_half = diag(Coef2)+ diag(Coef3,1)+ diag(Coef1,-1); %%%tridiagonal
    
    
       %Vort_half = AA_half\K_half;
       Vort_half =TDMA(AA_half,K_half);
       VH = [Vort(1,j),Vort_half',Vort(end,j)];  %%% ΠΡΣΘΗΚΗ TOP & BOTTOM WALL
       VORT_HALF(:,j-1) = VH';
       Vort_Half = [ Vort(:,1),VORT_HALF,Vort(:,end)]; %%% TELIKO GIA HALF TIME 
       
    
end %%%% END OF HALF TIME VORT CALCULATIONS %%%%
Vort = Vort_Half;    %%%% Redefine Vorticity from half time calcuklatuons








%%%%%% FULL TIME n +1 MARCHING IN X-DIRECTION %%%%%%

    % FOR VORTICITY
    Vort(:,1) = 2*(Psi(:,1)-Psi(:,2))/(StepX^2);                    %%% LEFT WALL
    Vort(:,end) = 2*(Psi(:,end)-Psi(:,end-1))/(StepX^2);            %%% RIGHT WALL
    Vort(end,:) = 2*(Psi(end,:)-Psi(end-1,:))/(StepY^2);            %%% BOTTOM WALL
    Vort(1,:) = (2*(Psi(1,:)-Psi(2,:))/(StepY^2)) - (2*U/StepY);    %%% TOP WALL


    %%% Vort matrix coefficients constants %%%
    C = (dt/2)*(1/Re)*(1/(StepY^2));
    D = (dt/2)*(1/Re)*(1/(StepX^2));
    
for i = 2:Ypoints-1 %%% ΓΡΑΜΜΕΣ
    
     
for j = 2:Xpoints -1 %%% ΣΤΗΛΕΣ

       %%% Variable Coeficients %%% (Τo i-1 γιατι το j ξεκινα απο 2!!)
       A(j-1) = (dt/2)*(1/(StepX*StepY))*(Psi(i+1,j) - Psi(i,j)); %%% HERE A & B IS j-1!!! we march in X
       B(j-1) = (dt/2)*(1/(StepX*StepY))*(Psi(i,j+1) - Psi(i,j)); 
    
       
   if j == 2
   K_full(1,j-1) = C*Vort(i-1,j) + (1-B(j-1)-2*C)*Vort(i,j) + (B(j-1)+C)*Vort(i+1,j) + (D)*Vort(i,j-1);    
   elseif j == Xpoints-1   
   K_full(1,j-1) = C*Vort(i-1,j) + (1-B(j-1)-2*C)*Vort(i,j) + (B(j-1)+C)*Vort(i+1,j) + (+D-A(j-1))*Vort(i,j+1);    
   else
   K_full(1,j-1) = C*Vort(i-1,j) + (1-B(j-1)-2*C)*Vort(i,j) + (B(j-1)+C)*Vort(i+1,j); 
   end
       
end
    %%% Coefficients of triadognal for full time %%%
    Coef4 = -D*ones(Ypoints - 3 , 1);
    Coef5 = (1-A+2*D)'.*ones(Ypoints - 2 , 1);
    Coef6 = (A(1:end-1)-D)'.*ones(Ypoints - 3 , 1);

     AA_full = diag(Coef5)+ diag(Coef6,1)+ diag(Coef4,-1); %%%tridiagonal
    
     %Vort_full = AA_full\K_full';
     Vort_full = TDMA( AA_full,K_full);
     VF = [Vort(i,1),Vort_full',Vort(i,end)];
     VORT_FULL(i-1,:) = VF';     % i-1 due to i starts from 2!!!
     Vort_Full = [  Vort(1,:);VORT_FULL;Vort(end,:)];

end %%% END OF FULL VORTICITY CALCULATIONS
Vort = Vort_Full;    %%%% Redefine Vorticity from full time calculatuons

    
    %%%%% FOR PSI CALCULATIONS %%%%
    Psi_A = (StepY/StepX)^2;
    Psi_B = StepY^2;
    
    
    for j = 2:Xpoints-1
        
        
        for i = 2:Ypoints-1
 
%%% PSI EQUATION IMPLICIT SCHEME%%%            
%          if i == 2
%          K_half_Psi(i-1) = -Vort(i,j)*Psi_B - Psi_A*(Psi(i,j+1) - 2*Psi(i,j) + Psi(i,j-1)) - Psi(i-1,j);
%          elseif i == Ypoints-1
%          K_half_Psi(i-1) = -Vort(i,j)*Psi_B - Psi_A*(Psi(i,j+1) - 2*Psi(i,j) + Psi(i,j-1)) - Psi(i+1,j);
%          else
%          K_half_Psi(i-1) = -Vort(i,j)*Psi_B - Psi_A*(Psi(i,j+1) - 2*Psi(i,j) + Psi(i,j-1));
%          end
%         
%         end
%         %%%% PSI half time coeficients %%%%
%        Coef1_Psi = ones(Ypoints - 3 , 1);  %%%%% LEFT 
%        Coef2_Psi = (-2)*ones(Ypoints - 2 , 1);  %%%%% CENTRAL
%        Coef3_Psi = ones(Ypoints - 3 , 1);  %%%%% RIGHT
%        
%        AA_half_Psi = diag(Coef2_Psi)+ diag(Coef3_Psi,1)+ diag(Coef1_Psi,-1); %%%tridiagonal
%         
%        Psi_half = TDMA(AA_half_Psi,K_half_Psi);
%        %Psi_half = AA_half_Psi\K_half_Psi;
%        PH = [Psi(1,j),Psi_half',Psi(end,j)];  %%% ΠΡΣΘΗΚΗ TOP & BOTTOM WALL
%        PSI_HALF(:,j-1) = PH';
%        Psi_Half = [ Psi(:,1), PSI_HALF,Psi(:,end)]; %%% TELIKO GIA HALF TIME 
       %%% PSI EQUATION EXPLICIT SCHEME%%%
            Psi(i,j) = (1/4)*(Vort(i,j)*(StepX^2) + (Psi(i,j+1)+Psi(i,j-1)) +  (Psi(i+1,j)+Psi(i-1,j)));
        end         
    end
    %Psi = Psi_Half;
 
    
    %%%% RELAXATION %%%%
    Vort = Vort_previous + Relax*(Vort-Vort_previous);
    Psi  = Psi_previous  + Relax*(Psi-Psi_previous);
    
    
    
     %%% Check for Converge %%%
    if time>dt*10
    ERROR(timer) = max(max(abs(Vort - Vort_previous)));
    if ERROR(timer) < 1e-4
        fprintf('man μου οι επαναληψεις είναι %d με βήμα %d\n', time/dt, dt);
        break
    elseif ERROR(timer) > 10
        
        disp('σκαει ο κωδικασ')
        break
    end
    timer=timer+1;
    end

end %%% END OF TIME %%%
toc

Array_vort(k,g) = Vort(Xpoints/2,Xpoints/2);
Array_Psi(k,g) = Psi(Xpoints/2,Xpoints/2);
end
end
%%%% PLOT GRID DEPENDENCE ANALYSIS %%%%
%%% Vorticity values %%%
figure
plot([30,40,50,100],Array_vort(:,1),'r','LineWidth',2)
hold on
plot([30,40,50,100],Array_vort(:,2),'b','LineWidth',2)

xlabel('Xpoints of Square Grid (Xpoints = Ypoints)')
ylabel('Vorticity values')
title('Implicit scheme Grid dependence (Στροβιλότητα)')
legend('Re = 1','Re = 100')
grid on

%%% Vorticity values %%%
figure
plot([30,40,50,100],Array_Psi(:,1),'r','LineWidth',2)
hold on
plot([30,40,50,100],Array_Psi(:,2),'b','LineWidth',2)

xlabel('Xpoints of Square Grid (Xpoints = Ypoints)')
ylabel('Psi function values')
title('Implicit scheme Grid dependence (Ροική Συνάρτηση)')
legend('Re = 1','Re = 100')
grid on
