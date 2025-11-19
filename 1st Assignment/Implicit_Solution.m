clc
clear
tic
lamdaX = 10; %% μηκος πλάκας
U = 2;       %%boundary condjtjon for x component jn m/s
lamdaY = 0.1; %% ύψος παρατήρησης

    
        
Xpoints = 4000;              %%%grjd pojnts jn x-djrectjon είναι οι στήλες j
Ypoints =400;             %%%grjd pojnts jn y-djrectjon είναι οι σειρές   i
StepX = (lamdaX/(Xpoints-1));   %%% step jn x-djrectjon
StepY = (lamdaY/(Ypoints-1));   %%% step jn y-djrectjon
Grid = zeros((Ypoints),(Xpoints)); %%%2D grjd

V = 1.4*10^(-5);        %%% Air's Viscocity
v = zeros(size(Grid));  %%% initial condition for v component
u = zeros(size(Grid));  %%% initial condition for u component  (έμπεριέχει No slip Condition)

u(:,1) = U;              %%%'Free stream'
u(end,:) = U;            %%% Free Stream
 
 
K = zeros(Ypoints-2,1);
Aa = zeros(Ypoints-2,1);
Bb = zeros(Ypoints-3,1);
 
%%%%% CRANK NICHOLSON METHOD  %%%%%

for j = 1:Xpoints-1
    
    for i = 2:Ypoints -1   
     
        
        %%% Constant Coefficients of the tridiagonal Matrix

      
        A = (1/StepX)*u(i,j);
        Aa(i-1) = A;                       %%%% Store here for Tridiagonal matrix
        B = ((0.25/StepY))*v(i,j);
        Bb(i-1) = B;                       %%%% Store here for Tridiagonal matrix
        
        C = V/(2*StepY^2);
        D = (A+2*C); %%% coeficient for u --> i
        E = (B-C);   %%% coeficient for u --> i+1
        Z = E;    %%% coeficient for u --> i-1
        
        if i==2
            K(i-1,1) = (A-2*C)*u(i,j) +  (C-B)*u(i+1,j) + (B+C)*u(i-1,j) - E*u(i-1,j+1);      %%%% No slip Condition bondary
        elseif i == Ypoints-1
            K(i-1,1) =  (A-2*C)*u(i,j) + (C-B)*u(i+1,j) + (B+C)*u(i-1,j) - E*u(i+1,j+1);     %%% Free stream Boundary
        else 
            K(i-1,1) =  (A-2*C)*u(i,j) + (C-B)*u(i+1,j) + (B+C)*u(i-1,j);
            
        end
        
    end  
    
        %%% Matrix construction %%%
        CC = C*ones(Ypoints - 2 , 1);
        DD = (Aa + 2*CC);               %%% Στον τριδιαγώνιο πίνακα αλλάζει ο όρος ανά σειρά!!! (Μου έβγαλε την παναγία να το βρώ!!!)
        
        CC2 = C*ones(Ypoints - 3 , 1);
        ZZ = (Bb(2:end)-CC2);
        EE = (Bb(1:end-1)-CC2);

        AA=diag(DD)+ diag(EE,1)+ diag(ZZ,-1); %%%tridiagonal
        
        
        
       Uf = AA\K;
       UF = [0,Uf',2];
       u(:,j+1) = UF';
       
     
      
        %%%% 1st Equatiom PDE %%%%   %% theres mistake here
        for vi = 2:Ypoints-1
        v(vi,j+1) = (StepY/StepX)*(u(vi,j)-u(vi,j+1))+ v(vi+1,j+1);
        end
       
   
    
end

toc



cmin = 0; % Minimum value
cmax = 2; % Maximum value

%Create a heatmap using imagesc with the Jet colormap and custom color limits
colormap(jet);
imagesc(flip(u), [cmin, cmax]);
% Set axis labels and title
xlabel('Χ points of Grid');
ylabel('Y points of Grid');
title('U compomtent of Boundary layer in m/s');
colorbar





u_array = u;

%%%% COMPARE WITH ANALYTICAL BLASIUS SOLUTION PDE %%%%
b=1;
k= StepX:StepX:lamdaX+StepX;
delta = zeros(length(k),1);
delta1 = zeros(length(k),1);
delta2 = zeros(length(k),1);
Cf_blasius = zeros(length(k),1);
for  k = StepX:StepX:lamdaX+StepX
     
    Rex = U*k/V; %Reynold Number
    
    %%% Thickness %%%
    delta(b,1) = 5*k/(Rex^(1/2));
    %%% Displacemnt Thickness %%%
    delta1(b) = 1.721*k/(Rex^(1/2));
    %%% momentum thicknes Thickness %%%
    delta2(b) = 0.664*k/(Rex^(1/2));
     %%% Friction Coeff %%%
    Cf_blasius(b) = 0.664/(Rex^(1/2));
    
    b=b+1; %counter


end



targetValue = 1.975;
array=u;
% Calculate absolute differences
differences = abs(array - targetValue);

% Find the row indices where the absolute difference is minimized
[minDifferences, col_indices] = min(differences, [], 1);

% Find the column indices where the absolute difference is minimized
%Rows = find(differences == minDifferences);
Columns=col_indices;
Rows=0:StepX:lamdaX;

% Διάγραμμα για πάχος οριακού στρώματος αναλυτικά και αριθμητικά
figure;
plot(Rows/lamdaX, StepY*Columns/max(delta), 'g','LineWidth',2.2);
xlabel('Μήκος επίπεδης πλάκας σε x/L');
ylabel('Πάχος οριακού στρώματος σε y/δmax');
title('Οριακό στρώμα επίπεδης πλάκας  IMPLICIT SOLUTION');
grid on
hold on
plot(Rows/lamdaX,delta/max(delta),'b--','LineWidth',2.5)
legend('Οριακό στρώμα αριθμητική λύση', 'Οριακό στρώμα Blasius');
xlim([0 1]);
%Απόκλιση 
Div1=(abs((StepY*Columns(end)-delta(end))/delta(end)));
DivAbs1=(norm(StepY*Columns'-delta))/norm(delta)*100;

% Διάγραμμα για πάχος μετατόπισης οριακού στρώματος αναλυτικά και αριθμητικά
figure
x = 0:StepX:lamdaX; %%% ΓΙΑ ΤΑ ΜΕΤΡΑ ΤΙΣ ΠΛΑΚΑΣ
disp =sum(1 - u(:,1:1:Xpoints)/U)*StepY; %% ολοκλήρωμα μετατόπισης

plot(x/lamdaX,(sum(1 - u(:,1:1:Xpoints)/U)*StepY)/max(delta),'c','LineWidth',2.2);
hold on
plot(x/lamdaX , delta1/max(delta),'k--','LineWidth',2.2)
xlabel('Μήκος επίπεδης πλάκας σε x/L');
ylabel('Πάχος μετατόπισης οριακού στρώματος μετατόπισης σε δ*/δmax');
title('δ* οριακού στρώματος IMPLICIT SOLUTION');
legend('Οριακό στρώμα αριθμητική λύση', 'Οριακό στρώμα Blasius')
grid on
%Απόκλιση
Div2 = abs((disp(end)-delta1(end))/delta1(end));
DivAbs2=(norm(StepY*Columns'-delta1))/norm(delta1)*100;

% Διάγραμμα για πάχος ΟΡΜΗΣ οριακού στρώματος αναλυτικά και αριθμητικά
figure
Theta = sum(u/U - (u.^2)/(U^2))*StepY;  %Ολοκλήρωμα ορμής
plot(x/lamdaX,(sum(u/U - (u.^2)/(U^2))*StepY)/max(delta),'c','LineWidth',2.2)
hold on
plot(x/lamdaX , delta2/max(delta),'m--','LineWidth',2.2)
xlabel('Μήκος επίπεδης πλάκας σε x/L');
ylabel('Πάχος ορμής οριακού στρώματος μετατόπισης σε δ2/δmax');
title('δ2 οριακού στρώματος IMPLICIT SOLUTION');
legend('Οριακό στρώμα αριθμητική λύση', 'Οριακό στρώμα Blasius')
grid on

                   %Απόκλιση
Div3 = abs((Theta(end)-delta2(end))/delta2(end));
DivAbs3=(norm(StepY*Columns'-delta2))/norm(delta2)*100;

% Διάγραμμα για τοπικό συντελεστή τριβής
Tw=2*Theta;                            %Διατμητική τάση

figure
plot(x/lamdaX,Cf_blasius,'r--','LineWidth',2.2)
hold on
plot(x/lamdaX,Tw./x,'b-','LineWidth',2)
xlabel('Μήκος επίπεδης πλάκας σε x/L');
ylabel('Cf');
title('Τοπικός συντελεστής τριβής IMPLICIT SOLUTION');
legend('Οριακό στρώμα αριθμητική λύση', 'Οριακό στρώμα Blasius')
grid on




fprintf('Η απόκλιση πάχους οριακού στρώματος για L=10m είναι%d %%\n', Div1*100);
fprintf('Η απόκλιση οριακού στρώματος μετατόπισης  για L=10m είναι %d %%\n', Div2*100);
fprintf('Η απόκλιση οριακού στρώματος ορμής  για L=10m είναι %d %%\n', Div3*100);


fprintf('Aπόλυτο σφάλμα πάχους οριακού στρώματος είναι%d %%\n', DivAbs1);
fprintf('Aπόλυτο σφάλμα οριακού στρώματος μετατόπισης είναι %d %%\n', DivAbs2);
fprintf('Aπόλυτο σφάλμα οριακού στρώματος ορμής είναι %d %%\n', DivAbs3);


