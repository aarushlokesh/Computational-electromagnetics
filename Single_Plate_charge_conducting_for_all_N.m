clc
clear all
% Define basic parameters
Z = 0;
k=1;
a  = input('Enter the length of the Plates');
Eo  = 8.854e-12;
for N=(1:1:60).^2 %total number of subsections
    for AA=1:length(N) 
    M  =  sqrt(N);
    b  = a/M; %Length of each subsection

                  

% Calculating source coordinates
P =  0; 
        
for p1 = 1:M
    for p2 = 1:M
        P  = P + 1;
        X(P) = b*(p1-0.5); %Calculate the co-ordinate along X direction
        Y(P) = b*(p2-0.5); %Calculate the co-ordinate along Y direction
    end
end



%Calculating the L matrix using the approximate method
 for i = 1:N
    for j = 1:N
        if(i==j)
            L(i,j) = ((b.*0.8814)/(pi.*Eo));
        else
            ll(i,j) = sqrt((X(j)-X(i))^2 + (Y(j)-Y(i))^2);
            L(i,j) = (((b./2).^2)/(pi.*Eo.*ll(i,j)));
        end
    end
 end

 % Potential Matrix 
for P = 1:N
      V(P) = 1; %Potential of each subsection
      
end
 
%Invert L and calculate sigma
 J = inv(L);
 sigma = J*V'; %Calculate the surface charge density
 

%Total Capacitance 
 Temp=sum(sigma(1:N));

C = Temp*(b.^2); %Calculate the total capacitance
c(k)=C; % Capacitances for various values of N
k=k+1;

    end
end

%plot charge density along subareas near centerline
midpt = ceil((N+1)/2);
if mod(N,2)==1  %N=odd case
    midpt = midpt - floor(M/2);
end
sigma3 = sigma(midpt:midpt+M-1);
xpts = 1/(2*M):1/M:1-1/(2*M);
figure, plot(xpts,sigma3,'o-')
xlabel('length');
ylabel('Charge density')

%Plot of Capacitance vs number of subsections
NN=1:1:length(c)
figure, plot(NN,c);
xlabel('Number of subsections in X and Y direction');
ylabel('Capaciatnce')



