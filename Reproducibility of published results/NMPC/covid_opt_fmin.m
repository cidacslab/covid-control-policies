
function cost=covid_opt(x,t,Model_0,N0,Usig,par)
global psir costd

% xv = x*ones(1,length(Usig));
% [minValue,closestIndex] = min(abs(Usig-xv));
% 
% % Calcula o valor 
% x = Usig(closestIndex);

[t,Psi,S,E,Ia,Is,H,U,R,D,Nw]=solver_model(x,t,Model_0,N0,par);

% Case 1
G1=5e5; %5e5
G2=1e5; %1e5 
G3=1e4; %1e4

% Case 2
G1=1e2; %5e5
G2=1e1; %1e5 
G3=1e3; %1e4


if t < 51
   nlhos = 466;
   nluti = 422;


for i=1:1:length(H)

if (H(i)>nlhos)||(U(i)>nluti)
    delta1(i) = Is'*Is + H'*H + U'*U;
else
    delta1(i)=0;
end
end

cost = Is(end) + G1*x + sum(delta1);

elseif t<193+21
    
nlhos = 1610;
nluti = 1210;

for i=1:1:length(H)

if (H(i)>nlhos)||(U(i)>nluti)
     delta2(i) = Is'*Is + H'*H + U'*U;
else
    delta2(i)=0;
end
end

       
        cost = Is(end) + G2*x + sum(delta2);

else

    nlhos = 1610;
    nluti = 1210;%-80;
    
    for i=1:1:length(H)

        if (H(i)>nlhos)||(U(i)>nluti)
             delta2(i) = Is'*Is + H'*H + U'*U;
        else
            delta2(i)=0;
        end
        end

        %Resultados 02/10: cost = Ia(end) + 20000*x + sum(delta2);
        cost = Is(end) + G3*x + sum(delta2);

end

if cost < costd
    psir = x;
    costd=cost;
end

end