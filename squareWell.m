function L = squareWell(a,N) 
% The program plots the wavefunctions for N[i,j], where i, j 
% are quantum numbers of a particle in a infinite well.
% The program assumes equal c_n's for each psi_n.
l = length(N);
E = zeros(l,1);
x = [-a/2:0.01:a/2]; % a well of length a can range from -a/2 to +a/2
lx = length(x);
p = zeros(l,lx); % l-row matrix. The row number is the quantum number n
P = ones(lx,1); % initializes the superposition wavefunction
pP = ones(lx,1); % initializes the PDF of P. pP = |<P|P>|^2
figure % draw in a new window
t = 0;
h = plot(x,pP,'LineWidth',2);
hold;
xlim([-a/2 a/2]) % correctly size the xaxis for the box.
%p1 = plot(x,p(1,:));
%p2 = plot(x,p(2,:));
for i = 1:l
    plt(i) = plot(x,p(i,:),'LineWidth',2);
end
legend('<\psi_a+\psi_b|\psi_a+\psi_b>','\psi_a','\psi_b')
plot(x,x*0,'Color','black','LineWidth',0.7);
set(gca,'FontSize',15)
while t >= 0
for i = 1:l
    E(i) = pi^2*N(i)^2/2/a^2;
    if mod(N(i),2) ~= 0
        p(i,:) = sqrt(2/a)*cos(pi*x*N(i)/a) * exp(-1i*E(i)*t);
    else
        p(i,:) = sqrt(2/a)*sin(pi*x*N(i)/a) * exp(-1i*E(i)*t);
    end
    %P = P.*p(i,:)';
    set(plt(i),'Ydata',real(p(i,:)));
    %set(p2,'Ydata',real(p(2,:)));
end
P = 1/sqrt(2)*(p(1,:)+p(2,:))';
pP = P.*conj(P); 
ylim([-2 4]) % uncomment this and last line and comment the next line and 
% l.42 to plot <x+y|x+y> instead of real(<y|x>)
%pP = p(1,:).*conj(p(2,:));
set(h,'Ydata',real(pP));
drawnow
%ylim([-3.5 3.5]);
%S_Psi = sum(P.*conj(P))*0.01 % Uncomment to verify <x+y|x+y> = 1
pause(.1)
t = t+0.01;
end
end
