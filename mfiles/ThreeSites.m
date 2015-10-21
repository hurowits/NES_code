clear all
tic
nu_vec=linspace(0,10,1000);
for iNu = 1:length(nu_vec)
    nu = nu_vec(iNu);
    W = zeros(3);
    W(1,2) = 1;
    W(2,3) = 1 ;
    W(3,1) = 1+nu;
    
    W(2,1) = 0;
    W(3,2) = 0;
    W(1,3) = 0;
    
    
     
    W(1,1) = -W(3,1) - W(2,1);
    W(2,2) = -W(3,2) - W(1,2);
    W(3,3) = -W(1,3) - W(2,3);
    
    A = diag(diag(W));
    B = W-A;
    
    A = 0.5*(W+W');
    B = 0.5*(W-W');
    lambda(iNu,:) = eig(W);
    gamma(iNu) = 2*trace(A^2+B^2)/trace(A)^2;
    detW(iNu)=det(W);
    detW2(iNu) = (trace(W)^3-3*trace(W)*trace(W^2)+2*trace(W^3))/6;
    detW3(iNu) = det(diag(diag(W)));
    
    norm_B(iNu) = B(1,2)^2+B(1,3)^2+B(2,3)^2;
end
%%
figure;
plot(gamma, imag(lambda))
xlabel('2var(\gamma)/<\gamma>^2');
ylabel('Im[\lambda]')

figure;
plot(gamma, real(lambda))
xlabel('2var(\gamma)/<\gamma>^2');
ylabel('Re[\lambda]')

%%
figure;
plot(nu_vec, imag(lambda))
xlabel('\nu');
ylabel('Im[\lambda]')

figure;
plot(nu_vec, real(lambda))

xlabel('\nu')
ylabel('Re[\lambda]')

%%

figure;
plot(real(lambda), imag(lambda),'.','MarkerSize',10)
xlabel('Re[\lambda]');
ylabel('Im[\lambda]');
%%
figure;
plot(nu_vec,gamma);
xlabel('\nu');
ylabel('\gamma');
%%
%%
figure;
plot(norm_B,real(lambda),norm_B,imag(lambda));
xlabel('|B|');
ylabel('\lambda');

toc