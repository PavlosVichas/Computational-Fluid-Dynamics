clc
clear
tic
lamdaX = 10; %% μηκος πλάκας
U = 2;       %%boundary condjtjon for x component jn m/s
lamdaY = 0.1; %% ύψος παρατήρησης


Xpoints = 5000 ;              %%%grjd pojnts jn x-djrectjon  είναι στήλες
Ypoints = 150  ;             %%%grjd pojnts jn y-djrectjon  είναι γραμμές
StepX = (lamdaX/(Xpoints-1));   %%% step jn x-djrectjon
StepY = (lamdaY/(Ypoints-1));   %%% step jn y-djrectjon
Grid = zeros((Ypoints),(Xpoints)); %%%2D grjd

V = 1.4*10^(-5);        %%% Air's Viscocity
v = zeros(size(Grid));  %%% initial condition for v component
u = zeros(size(Grid));  %%% initial condition for u component 
              
u(end,:) = 0;           %%% No slip condition 
u(:,1) = U;             %%% Free stream στήλη
u(1,:)=U;               %%% Free stream στο άπειρο

%%%% i for Y direction %%%%   %%% i είναι οι σειρές %%%
%%%% j for X direction %%%%   %%% j ειναι οι στήλες %%%

for  j = 1:Xpoints-1             
    
   for i = Ypoints-1:-1:2
    
    %%% Second   equation
    du2dy2 = V*((u(i-1,j) - 2*(u(i,j)) + u(i+1,j))/(StepY^2));
    vdudy = v(i,j)*(u(i-1,j) - u(i,j))/StepY;
    u(i,j+1) = (du2dy2 - vdudy)*(StepX/u(i,j)) + u(i,j);
    
    %%% First equation
    dudx = (u(i,j+1)-u(i,j))/StepX;
    v(i-1,j) = (-dudx)*StepY + v(i,j);
    
   end
 
end

u_array = flip(u);

figure;
% Define custom color limits
cmin = 0; % Minimum value
cmax = 2; % Maximum value
%Create a heatmap using imagesc with the Jet colormap and custom color limits
colormap(jet);
imagesc(u, [cmin, cmax]);
% Set axis labels and title
xlabel('Χ points of Grid');
ylabel('Y points of Grid');
title('U compomtent of Boundary layer in m/s');
colorbar
toc




%%%% COMPARE WITH ANALYTICAL BLASIUS SOLUTION PDE %%%%
b=1;
k=0:0.01:lamdaX;
delta = zeros(length(k),1);
delta1 = zeros(length(k),1);
delta2 = zeros(length(k),1);
Cf_blasius = zeros(length(k),1);

for  k = 0:0.01:lamdaX
     
    Rex = U*k/V;
    
    %%% Thickness %%%
    delta(b,1) = 5*k/(Rex^(1/2));
    %%% Displacemnt Thickness %%%
    delta1(b) = 1.721*k/(Rex^(1/2));
    %%% momentum thicknes Thickness %%%
    delta2(b) = 0.664*k/(Rex^(1/2));
    %%% Friction Coeff %%%
    Cf_blasius(b) = 0.664/(Rex^(1/2));
    
    b=b+1;


end
k=0:0.01:lamdaX;

%%% ΕΥΡΕΣΗ ΠΑΧΟΥΣ ΟΡΙΑΚΟΥ ΣΤΡΩΜΑΤΟΣ %%%
u_array = u_array';
[rows, cols] = size(u_array);
% Find  values between 1.98 and 1.99
values = find(u_array >= 1.9815 & u_array <= 1.9830);
% Extract row and column indices for each valid value
Rows = mod(values - 1, rows) + 1; % Convert index to row number
Columns = floor((values - 1) / rows) + 1; % Convert index to column number

% Διάγραμμα για πάχος οριακού στρώματος αναλυτικά και αριθμητικά
figure;
plot((Rows*StepX)/lamdaX, StepY*Columns/max(delta), 'g','LineWidth',2.2);
xlabel('Μήκος επίπεδης πλάκας σε x/L');
ylabel('Πάχος οριακού στρώματος σε y/δ max');
title('Οριακό στρώμα επίπεδης πλάκας,EXPLICIT SOLUTION');
grid on
hold on
plot(k/lamdaX,delta/max(delta),'b--','LineWidth',2.5)
legend('Οριακό στρώμα με Αριθμιτική λύση', 'Οριακό στρώμα Blasius');
xlim([0 1])
                %Απόκλιση 
Div1=abs((StepY*Columns(end)-delta(end))/delta(end));

% Διάγραμμα για πάχος μετατόπισης οριακού στρώματος αναλυτικά και αριθμητικά
figure
x = 0:0.002:10-0.002; 
disp =sum(1 - u(:,1:1:5000)/U)*StepY; %% ολοκλήρωμα μετατόπισης

plot(x/lamdaX,(sum(1 - u(:,1:1:5000)/U)*StepY)/max(delta),'c','LineWidth',2.2);
hold on
plot(k/lamdaX , delta1/max(delta),'k--','LineWidth',2.2)
xlabel('Μήκος επίπεδης πλάκας σε x/L');
ylabel('Πάχος οριακού στρώματος μετατόπισης σε δ*/max δ');
title('δ* οριακού στρώματος EXPLICIT SOLUTION');
legend('Οριακό στρώμα με Αριθμιτική λύση', 'Οριακό στρώμα Blasius')
grid on
                %Απόκλιση
Div2 = abs((disp(end)-delta1(end))/delta1(end));

% Διάγραμμα για πάχος ΟΡΜΗΣ οριακού στρώματος αναλυτικά και αριθμητικά
figure
Theta = sum(u/U - (u.^2)/(U^2))*StepY;  %Ολοκλήρωμα ορμής 
plot(x/lamdaX,(sum(u/U - (u.^2)/(U^2))*StepY)/max(delta),'c','LineWidth',2.2)
hold on
plot(k/lamdaX , delta2/max(delta),'m--','LineWidth',2.2)
xlabel('Μήκος επίπεδης πλάκας σε x/L');
ylabel('Πάχος ορμής οριακού στρώματος σε δ2/δmax');
title('δ2 οριακού στρώματος EXPLICIT SOLUTION');
legend('Οριακό στρώμα με Αριθμιτική λύση', 'Οριακό στρώμα Blasius')
grid on
                   %Απόκλιση
Div3 = abs((Theta(end)-delta2(end))/delta2(end));

% Διάγραμμα για τοπικό συντελεστή τριβής
Tw=2*Theta; %Διατμητική τάση

figure
plot(k,Cf_blasius,'r--','LineWidth',2.2)
hold on
plot(x,Tw./x,'b-','LineWidth',2)
xlabel('Μήκος επίπεδης πλάκας σε m');
ylabel('Cf');
title('Τοπικός συντελεστής τριβής, EXPLICIT SOLUTION');
legend('Οριακό στρώμα με Αριθμιτική λύση', 'Οριακό στρώμα Blasius')
grid on




fprintf('Η απόκλιση πάχους οριακού στρώματος είναι%d %%\n', Div1*100);
fprintf('Η απόκλιση οριακού στρώματος μετατόπισης είναι %d %%\n', Div2*100);
fprintf('Η απόκλιση οριακού στρώματος ορμής είναι %d %%\n', Div3*100);