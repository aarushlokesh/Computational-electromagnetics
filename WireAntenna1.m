clc;
clear all; 
%Parameters
c = 3e8; % speed of light
f = 300e6; %frequency
lambda = c/f; % wavelength
W = 2*pi*(f); % angular frequency
L = 0:0.01:2*lambda; % length of wire
N = input('Enter number of subsections:');
k = 2*pi/lambda; % wavenumber
mu_0 = (4*pi)*1e-7; % permeability of free space
eps_0 = 1e-9/(36*pi); % permittivity of free space

for i = 1:length(L)
    l = L(i);
    delta = l/N;
    a = l/(2*74.2);
    
    %Calculate Psi
    xi = zeros(N); %Psi(n,m)
    xi1 = zeros(N); %Psi(n-plus, m-plus)
    xi2 = zeros(N); %Psi(n-minus, m-plus)
    xi3 = zeros(N); %Psi(n-plus, m-minus)
    xi4 = zeros(N); %Psi(n-minus, m-minus)
    
    for n = 1:N
        for m = 1:N
            
            %Psi(n,m)
            if n == m
                Rmn = a;
                xi(n,m) = 1/(2*pi*delta)*log(delta/Rmn)-1i*k/(4*pi);
            elseif n < m
                Rmn = sqrt(a^2+(((m)*delta+delta/2)-((n)*delta+delta/2))^2);
                xi(n,m) = exp(-1i*k*Rmn)/(4*pi*Rmn);
            else
                Rmn = sqrt(a^2 +(((n-1)*delta+delta/2)-((m-1)*delta+delta/2))^2);
                xi(n,m) = exp(-1i*k*Rmn)/(4*pi*Rmn);
            end
            
            %Psi(n-plus, m-plus)
            if n == m
                Rmn = a;
                xi1(n,m) = 1/(2*pi*delta)*log(delta/Rmn)-1i*k/(4*pi);
            elseif n < m
                Rmn = sqrt(a^2+(m*delta-n*delta)^2);
                xi1(n,m) = exp(-1i*k*Rmn)/(4*pi*Rmn);
            else
                Rmn = sqrt(a^2+(n*delta-m*delta)^2);
                xi1(n,m) = exp(-1i*k*Rmn)/(4*pi*Rmn);
            end

            %Psi(n-minus, m-plus)
            if n == m
                Rmn = sqrt(a^2+delta^2);
                xi2(n,m) = 1/(2*pi*delta)*log(delta/Rmn)-1i*k/(4*pi);
            elseif n < m
                Rmn = sqrt(a^2+(m*delta-(n-1)*delta)^2);
                xi2(n,m) = exp(-1i*k*Rmn)/(4*pi*Rmn);
            else
                Rmn = sqrt(a^2 + (n*delta-m*delta)^2);
                xi2(n,m) = exp(-1i*k*Rmn)/(4*pi*Rmn);
            end

            %Psi(n-plus, m-minus)
            if n == m
                Rmn = sqrt(a^2+delta^2);
                xi3(n,m) = 1/(2*pi*delta)*log(delta/Rmn)-1i*k/(4*pi);
            elseif n < m
                Rmn = sqrt(a^2+(m*delta-n*delta)^2);
                xi3(n,m) = exp(-1i*k*Rmn)/(4*pi*Rmn);
            else 
                Rmn = sqrt(a^2 +(n*delta-(m-1)*delta)^2);
                xi3(n,m) = exp(-1i*k*Rmn)/(4*pi*Rmn);
            end

            %Psi(n-minus, m-minus)
            if n == m
                Rmn = a;
                xi4(n,m) = 1/(2*pi*delta)*log(delta/Rmn)-1i*k/(4*pi);
            elseif n < m
                Rmn = sqrt(a^2+((m-1)*delta-(n-1)*delta)^2);
                xi4(n,m) = exp(-1i*k*Rmn)/(4*pi*Rmn);
            else
                Rmn = sqrt(a^2+((n-1)*delta-(m-1)*delta)^2);
                xi4(n,m) = exp(-1i*k*Rmn)/(4*pi*Rmn);
            end  
        end
    end

    %Calculate Impedance Matrix
    Z = (1i*W*mu_0*delta)*(delta*xi)+((1/(1i*W*eps_0))*(xi1-xi2-xi3-xi4));
    Y = inv(Z);

    for j = 1:N
        if j == ceil(N/2)
            V(j) = 1;
        else
            V(j) = 0;
        end
    end

I = V/Z;
Yin(i) = Y(length(N/2), length(N/2));
end

Impedence = 1/max(Yin);

figure(1)
plot(L,real(Yin)/1e-3); 
grid on;
xlabel('length of antenna')
ylabel(' G in millimhos')

figure(2)
plot(L,imag(Yin)/1e-3); 
grid on;
xlabel('length of antenna')
ylabel(' B in millimhos')
