function [I0,I90,I45,I135] = shearing_interferometry(lambda,d1,L0,d2,Uo,D,delta_s,f,N0,delta_x)
%% Front surface of metalens %%
x0 = linspace(-L0/2, L0/2,N0);
y0 = linspace(-L0/2, L0/2,N0);
[X0, Y0] = meshgrid(x0, y0);

L1=N0*lambda*d1/L0;
k=2*pi/lambda;

x1 = linspace(-L1/2, L1/2,N0);
y1 = linspace(-L1/2, L1/2,N0);
[X1, Y1] = meshgrid(x1, y1);

F0=exp(1i*k*d1)/(1i*lambda*d1)*exp(1i*k/2/d1*(X1.^2+Y1.^2));
F=exp(1i*k/2/d1*(X0.^2+Y0.^2));
FU=(L0*L0/N0/N0).*fftshift(fft2(Uo.*F));
U1=F0.*FU;

%% Diaphragm %%
pxy=zeros(N0,N0);
for n=1:N0
    for m=1:N0
        if x1(n)^2+y1(m)^2<=(D/2)^2
            pxy(n,m)=1;
        end
    end
end
pxy=ones(N0,N0);

%% Images on CCD %%
L2  = N0*lambda*d2/L1;
x2 = linspace(-L2/2, L2/2,N0);
y2 = linspace(-L2/2, L2/2,N0);
[X2, Y2] = meshgrid(x2, y2);

E_back_lens_LCP1_L = U1.*pxy.*exp(-1i*k.*(sqrt(((X1-delta_s/2+delta_x).^2+Y1.^2)+f^2)-f));
E_back_lens_LCP2_L = U1.*pxy.*exp(-1i*k.*(sqrt(((X1-delta_s/2+delta_x).^2+Y1.^2)+f^2)-f)).*1i;
E_back_lens_RCP1_L = U1.*pxy.*exp(-1i*k.*(sqrt(((X1+delta_s/2+delta_x).^2+Y1.^2)+f^2)-f));
E_back_lens_RCP2_L = -U1.*pxy.*exp(-1i*k.*(sqrt(((X1+delta_s/2+delta_x).^2+Y1.^2)+f^2)-f)).*1i;

E_back_lens_LCP1_R = U1.*pxy.*exp(-1i*k.*(sqrt(((X1-delta_x).^2+(Y1-delta_s/2).^2)+f^2)-f));
E_back_lens_LCP2_R = U1.*pxy.*exp(-1i*k.*(sqrt(((X1-delta_x).^2+(Y1-delta_s/2).^2)+f^2)-f)).*1i;
E_back_lens_RCP1_R = U1.*pxy.*exp(-1i*k.*(sqrt(((X1-delta_x).^2+(Y1+delta_s/2).^2)+f^2)-f));
E_back_lens_RCP2_R = -U1.*pxy.*exp(-1i*k.*(sqrt(((X1-delta_x).^2+(Y1+delta_s/2).^2)+f^2)-f)).*1i;


F0=exp(1i*k*d2)/(1i*lambda*d2)*exp(1i*k/2/d2*(X2.^2+Y2.^2));
F=exp(1i*k/2/d2*(X1.^2+Y1.^2));
E_CCD_LCP1_L=(L1*L1/N0/N0).*fft2(E_back_lens_LCP1_L.*F).*F0;
E_CCD_LCP2_L=(L1*L1/N0/N0).*fft2(E_back_lens_LCP2_L.*F).*F0;
E_CCD_RCP1_L=(L1*L1/N0/N0).*fft2(E_back_lens_RCP1_L.*F).*F0;
E_CCD_RCP2_L=(L1*L1/N0/N0).*fft2(E_back_lens_RCP2_L.*F).*F0;

E_CCD_LCP1_R=(L1*L1/N0/N0).*fft2(E_back_lens_LCP1_R.*F).*F0;
E_CCD_LCP2_R=(L1*L1/N0/N0).*fft2(E_back_lens_LCP2_R.*F).*F0;
E_CCD_RCP1_R=(L1*L1/N0/N0).*fft2(E_back_lens_RCP1_R.*F).*F0;
E_CCD_RCP2_R=(L1*L1/N0/N0).*fft2(E_back_lens_RCP2_R.*F).*F0;

ExL = E_CCD_LCP1_L + E_CCD_RCP1_L;
ExR = E_CCD_LCP1_R+E_CCD_RCP1_R;
I0 = ExL.*conj(ExL)+ExR.*conj(ExR);
EyL = E_CCD_LCP2_L + E_CCD_RCP2_L;
EyR = E_CCD_LCP2_R+E_CCD_RCP2_R;
I90 = EyL.*conj(EyL)+EyR.*conj(EyR);
S_45 = [1 0; 0 0] * sqrt(2)/2 * [1 1;-1 1] ;
S_135 = [0 0; 0 1] * sqrt(2)/2 * [1 1;-1 1] ;

for indx = 1:length(x2)
    for indy = 1:length(y2)

        E_CCD_45L(indx, indy) = sum(S_45 * [ExL(indx,indy) ; EyL(indx,indy)]);
        E_CCD_45R(indx, indy) = sum(S_45 * [ExR(indx,indy) ; EyR(indx,indy)]);
        E_CCD_135L(indx, indy) = sum(S_135 * [ExL(indx,indy) ; EyL(indx,indy)]);
        E_CCD_135R(indx, indy) = sum(S_135 * [ExR(indx,indy) ; EyR(indx,indy)]);

    end
end

I45 = E_CCD_45L.*conj(E_CCD_45L)+E_CCD_45R.*conj(E_CCD_45R);
I135 = E_CCD_135L.*conj(E_CCD_135L)+E_CCD_135R.*conj(E_CCD_135R);

end