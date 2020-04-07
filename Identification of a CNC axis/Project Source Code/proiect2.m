close all
t=double(pop_denis.X.Data)';
d=double(pop_denis.Y(1,3).Data)'; %semnal de comanda ce preciseaza sensul
w=double(pop_denis.Y(1,2).Data)'; %viteza unghiulara
theta=double(pop_denis.Y(1,1).Data)'; %pozitia unghiulara

figure
subplot(3,1,1); plot(t,d); title('Intrarea (u)'); xlabel('t[s]'); ylabel('Factor de umplere PWM');
subplot(3,1,2); plot(t,w); title('Viteza unghiulara (w)'); xlabel('t[s]'); ylabel('w[rad/s]');
subplot(3,1,3); plot(t,theta); title('Pozitia (theta)'); xlabel('t[s]'); ylabel('pozitia[nr de impulsuri]');

Te=t(2)-t(1);
i1=1034;
i2=1740;
i3=4597;
i4=5295;
i5=8386;
i6=9080;
wi=w;
wi(i1:i2)=interp1([t(i1) t(i2)],[wi(i1) wi(i2)],t(i1:i2));
wi(i3:i4)=interp1([t(i3) t(i4)],[wi(i3) wi(i4)],t(i3:i4));
wi(i5:i6)=interp1([t(i5) t(i6)],[wi(i5) wi(i6)],t(i5:i6));
% 
% plot(t,w); title('Viteza unghiulara (w)'); xlabel('t[s]'); ylabel('w[rad/s]');
% plot(t,wi); title('Viteza unghiulara interpolata(wi)'); xlabel('t[s]'); ylabel('wi[rad/s]');
w=wi;

figure
subplot(3,1,1);plot(t,d);
subplot(3,1,2);plot(t,w);
subplot(3,1,3);plot(t,theta);

t1=1745;
t2=4596;
t3=5338;
t4=8256;
data_id_w=iddata(w(t1:t2),d(t1:t2),Te);%date identificare viteza
data_id_th=iddata(theta(t1:t2),w(t1:t2),Te); %date identificare pozitie
data_v_w=iddata(w(t3:t4),d(t3:t4),Te);%date validare viteza
data_v_th=iddata(theta(t3:t4),w(t3:t4),Te); %date validare pozitie
data_g_w=iddata(w,d,Te);%date generale la viteza
data_g_th=iddata(theta,w,Te); %date generale pozitie


%% ARX

model_arx_w=arx(data_id_w,[1 1 1]);
Hd_w_d_arx=tf(model_arx_w.B,model_arx_w.A,Te,'variable','z^-1');%functia de transfer de la intrare d cu iesire w
model_arx_th=arx(data_id_th,[1 1 1]);
Hd_th_w_arx=tf(model_arx_th.B,model_arx_th.A,Te,'variable','z^-1');%functia de transfer de la intrare w cu iesire th

H_w_d_arx=d2c(Hd_w_d_arx,'zoh');

figure
resid(model_arx_w,data_v_w,'corr',5);
figure
compare(model_arx_w,data_g_w)

figure
resid(model_arx_th,data_v_th,'corr',5);
figure
compare(model_arx_th,data_g_th)


%% ARMAX
%w
model_armax_w=armax(data_id_w,[1 1 1 1]);
Hd_w_d_armax=tf(model_armax_w.B,model_armax_w.A,Te,'variable','z^-1');%functia de transfer de la intrare d cu iesire w

figure
resid(model_armax_w,data_v_w,'corr',5);
figure
compare(model_armax_w,data_g_w)

%th
model_armax_th=armax(data_id_th,[1 1 1 1]);
Hd_th_w_armax=tf(model_armax_th.B,model_armax_th.A,Te,'variable','z^-1');%functia de transfer de la intrare d cu iesire w

figure
resid(model_armax_th,data_v_th,'corr',5);
figure
compare(model_armax_th,data_g_th)

H_w_d_arx=d2c(Hd_w_d_armax,'zoh');
H_w_d_armax=d2c(Hd_w_d_armax,'zoh');


H_th_w_arx=d2c(Hd_th_w_armax,'zoh');
H_th_w_armax=d2c(Hd_th_w_armax,'zoh');


%% IV4
%w
m_iv_w=iv4(data_id_w,[1 1 1]);
Hd_w_d_iv=tf(m_iv_w.B,m_iv_w.A,Te,'variable','z^-1');
H_w_d_iv=d2c(Hd_w_d_iv,'zoh')

figure
resid(m_iv_w,data_v_w,'corr',5);
figure
compare(m_iv_w,data_g_w)

%th
model_iv_th=iv4(data_id_th,[1 1 1]);
Hd_th_w_iv=tf(model_iv_th.B,model_iv_th.A,Te,'variable','z^-1');%functia de transfer de la intrare d cu iesire w

figure
resid(model_iv_th,data_v_th,'corr',5);
figure
compare(model_iv_th,data_g_th)


%% OE ->Good
%w
m_oe_w=oe(data_id_w,[1 1 1]);
Hd_w_d_oe=tf(m_oe_w.B,m_oe_w.F,Te,'variable','z^-1');
H_w_d_oe=d2c(Hd_w_d_oe,'zoh')
figure
resid(m_oe_w,data_v_w,'corr',30);
figure
compare(m_oe_w,data_g_w)

%th
model_oe_th=oe(data_id_th,[1 1 1]);
Hd_th_w_oe=tf(model_oe_th.B,model_oe_th.F,Te,'variable','z^-1');%functia de transfer de la intrare d cu iesire w

figure
resid(model_oe_th,data_v_th,'corr',30);
figure
compare(model_oe_th,data_g_th)
%% Modele in continuu iv+oe
H_w_d_iv=d2c(Hd_w_d_iv,'zoh')
H_th_w_iv=tf([5.1],[1 0])

H_w_d_oe=d2c(Hd_w_d_oe,'zoh')
H_th_w_oe=tf([5.0850],[1 0]) %5.0850=indice z^-1 numarator / Te



%% Ex 1------------------------------------------- Proiect in caz de nu obt un model valid (PEM)
%ARX cu PEM
m_arx_w_pem=pem(data_id_w,model_arx_w);
figure
resid(m_arx_w_pem,data_v_w,'corr',10);
figure
compare(m_arx_w_pem,data_g_w)

m_arx_th_pem=pem(data_id_th,model_arx_th);
figure
resid(m_arx_th_pem,data_v_th,'corr',10);
figure
compare(m_arx_th_pem,data_g_th)

%% ARMAX cu PEM
m_arx_w_pem=pem(data_id_w,model_arx_w);
figure
resid(m_arx_w_pem,data_v_w,'corr',5);
figure
compare(m_arx_w_pem,data_g_w)

m_arx_th_pem=pem(data_id_th,model_arx_th);
figure
resid(m_arx_th_pem,data_v_th,'corr',5);
figure
compare(m_arx_th_pem,data_g_th)


%% ex 2 OE cu PEM
m_armax_w_pem=pem(data_id_w,model_armax_w);
figure
resid(m_armax_w_pem,data_v_w,'corr',5);
figure
compare(m_armax_w_pem,data_g_w)

m_armax_th_pem=pem(data_id_th,model_armax_th);
figure
resid(m_armax_th_pem,data_v_th,'corr',5);
figure
compare(m_armax_th_pem,data_g_th)


%% Ex3 DECIMARE
i7=2906;
i8=2910;
N=i8-i7+1;
t_dec=t(1:N:end);
d_dec=d(1:N:end);
w_dec=w(1:N:end);
th_dec=theta(1:N:end);
plot(t_dec,w_dec);
Te_dec=Te*N;
t1d=round(t1/N);
t2d=round(t2/N);
t3d=round(t3/N);
t4d=round(t4/N);

data_id_w_dec=iddata(w_dec(t1d:t2d),d_dec(t1d:t2d),Te_dec);%date identificare
data_id_th_dec=iddata(th_dec(t1d:t2d),w_dec(t1d:t2d),Te_dec);
data_v_w_dec=iddata(w_dec(t3d:t4d),d_dec(t3d:t4d),Te_dec);%date validare
data_v_th_dec=iddata(th_dec(t3d:t4d),w_dec(t3d:t4d),Te_dec);
data_g_w_dec=iddata(w_dec,d_dec,Te_dec);%date generale
data_g_th_dec=iddata(th_dec,w_dec,Te_dec);

%% Ex3 Arx cu decimare
model_arx_w_dec=arx(data_id_w_dec,[1 1 1]);
Hd_w_d_arx_dec=tf(model_arx_w_dec.B,model_arx_w_dec.A,Te_dec,'variable','z^-1');%functia de transfer de la intrare d cu iesire w
model_arx_th_dec=arx(data_id_th_dec,[1 1 1]);
Hd_th_w_arx=tf(model_arx_th_dec.B,model_arx_th_dec.A,Te_dec,'variable','z^-1');%functia de transfer de la intrare w cu iesire th

H_testare=d2c(Hd_w_d_arx_dec,'zoh')

figure
resid(model_arx_w_dec,data_v_w_dec,'corr',20);
figure
compare(model_arx_w_dec,data_g_w_dec)

figure
resid(model_arx_th_dec,data_v_th_dec,'corr',20);
figure
compare(model_arx_th_dec,data_g_th_dec)

%cu PEM si decimare
m_arx_w_pem_dec=pem(data_id_w_dec,model_arx_w_dec);

figure
resid(m_arx_w_pem_dec,data_v_w_dec,'corr',5);
figure
compare(m_arx_w_pem_dec,data_g_w_dec)

m_arx_th_pem_dec=pem(data_id_th_dec,model_arx_th_dec);
figure
resid(m_arx_th_pem_dec,data_v_th_dec,'corr',5);
figure
compare(m_arx_th_pem_dec,data_g_th_dec)

%% ARMAX cu decimare
model_armax_w_dec=armax(data_id_w_dec,[1 1 1 1]);
Hd_w_d_arx_dec=tf(model_armax_w_dec.B,model_armax_w_dec.A,Te_dec,'variable','z^-1');%functia de transfer de la intrare d cu iesire w
model_armax_th_dec=armax(data_id_th_dec,[1 1 1 1]);
Hd_th_w_arx=tf(model_armax_th_dec.B,model_armax_th_dec.A,Te_dec,'variable','z^-1');%functia de transfer de la intrare w cu iesire th
Htest1=d2c(Hd_w_d_arx_dec,'zoh')
figure
resid(model_armax_w_dec,data_v_w_dec,'corr',10);
figure
compare(model_armax_w_dec,data_g_w_dec)
figure
resid(model_armax_th_dec,data_v_th_dec,'corr',10);
figure
compare(model_armax_th_dec,data_g_th_dec)

%cu PEM si decimare
m_armax_w_pem_dec=pem(data_id_w_dec,model_armax_w_dec);
Hd_w_d_armax_pem_dec=tf(m_armax_w_pem_dec.B,m_armax_w_pem_dec.A,Te_dec,'variable','z^-1');
H_w_d_armax_pem_dec=d2c(Hd_w_d_armax_pem_dec,'zoh')

figure
resid(m_armax_w_pem_dec,data_v_w_dec,'corr',20);
figure
compare(m_armax_w_pem_dec,data_g_w_dec)

m_armax_th_pem_dec=pem(data_id_th_dec,model_armax_th_dec);
figure
resid(m_armax_th_pem_dec,data_v_th_dec,'corr',20);
figure
compare(m_armax_th_pem_dec,data_g_th_dec)


%% Ex4 OE cu DECIMARE
m_oe_w_dec=oe(data_id_w_dec,[1 1 1]);
Hd_w_d_oe_dec=tf(m_oe_w_dec.B,m_oe_w_dec.F,Te_dec,'variable','z^-1');
figure
resid(m_oe_w_dec,data_v_w_dec,'corr',5);
figure
compare(m_oe_w_dec,data_g_w_dec)

model_oe_th_dec=oe(data_id_th_dec,[1 1 1]);
Hd_th_w_oe=tf(model_oe_th_dec.B,model_oe_th_dec.F,Te,'variable','z^-1');%functia de transfer de la intrare d cu iesire w

figure
resid(model_oe_th_dec,data_v_th_dec,'corr',5);
figure
compare(model_oe_th_dec,data_g_th_dec)

H_test1=d2c(Hd_w_d_oe_dec,'zoh')

%Decimare + PEM
m_oe_w_dec_pem_dec=pem(data_id_w_dec,m_oe_w_dec);
figure
resid(m_oe_w_dec_pem_dec,data_v_w_dec,'corr',5);
figure
compare(m_oe_w_dec_pem_dec,data_g_w_dec)

m_oe_th_pem_dec=pem(data_id_th_dec,model_oe_th_dec);
figure
resid(m_oe_th_pem_dec,data_v_th_dec,'corr',5);
figure
compare(m_oe_th_pem_dec,data_g_th_dec);


%% EX 5
data_id=iddata([theta(t1:t2),w(t1:t2)],d(t1:t2),Te);
data_v=iddata([theta(t3:t4),w(t3:t4)],d(t3:t4),Te);
data_g=iddata([theta,w],d,Te);
M_pem=pem(data_id);

figure
resid(M_pem,data_v,'corr',5);
figure
compare(M_pem,data_g);


%% Ex 6 d->theta direct
data_dw_id=iddata(theta(t1:t2),d(t1:t2),Te);
data_dw_v=iddata(theta(t3:t4),d(t3:t4),Te);
data_dw_g=iddata(theta,d,Te);
m_oe_dw=oe(data_dw_id,[2 2 1]);
Hd_dw_oe=tf(m_oe_dw.B,m_oe_dw.F,Te,'variable','z^-1');
Hd_dw=d2c(Hd_dw_oe,'zoh');

figure
resid(m_oe_dw,data_dw_v,'corr',5);
figure
compare(m_oe_dw,data_dw_g);

m_oe_dw_pem=pem(data_dw_id,m_oe_dw);
figure
resid(m_oe_dw_pem,data_dw_v,'corr',5);
figure
compare(m_oe_dw_pem,data_dw_g);

%% lab 12

%% #1
ssest(data_id_w,1:10);
nx_w=2;

ssest(data_id_th,1:10);
nx_th=1;

%% #2
m_n4sid_w=n4sid(data_id_w,nx_w);
m_n4sid_th=n4sid(data_id_th,nx_th);

figure
resid(m_n4sid_w,data_v_w,'corr',5);
figure
compare(m_n4sid_w,data_g_w);

figure
resid(m_n4sid_th,data_v_th,'corr',5);
figure
compare(m_n4sid_th,data_g_th);


m_ssest_w=ssest(data_id_w,1);
m_ssest_th=ssest(data_id_th,1);

figure
resid(m_ssest_w,data_v_w,'corr',5);
compare(m_ssest_w,data_g_w);

figure
resid(m_ssest_th,data_v_th,'corr',5);
compare(m_ssest_th,data_g_th);


%% #3 n4sid cu decimare
%ordin

ssest(data_id_w_dec,1:10);
nx_w=2;

ssest(data_id_th_dec,1:10);
nx_th=1;


%n4sid cu decimare
m_n4sid_w_dec=n4sid(data_id_w_dec,1);
m_n4sid_th_dec=n4sid(data_id_th_dec,1);

figure
resid(m_n4sid_w_dec,data_v_w_dec,'corr',5);
compare(m_n4sid_w_dec,data_g_w_dec);

figure
resid(m_n4sid_th_dec,data_v_th_dec,'corr',5);
compare(m_n4sid_th_dec,data_g_th_dec);

%ssest
m_ssest_w_dec=ssest(data_id_w_dec,1);
m_ssest_th_dec=ssest(data_id_th_dec,1);

figure
resid(m_ssest_w_dec,data_v_w,'corr',5);
compare(m_ssest_w_dec,data_g_w);

figure
resid(m_ssest_th_dec,data_v_th,'corr',5);
compare(m_ssest_th_dec,data_g_th);

%% #4 d->theta, +#1 +#2
    
ssest(data_id,1:10);
nx=2;

m_n4sid=n4sid(data_id,nx);

figure
resid(m_n4sid,data_v,'corr',5);
compare(m_n4sid,data_g);

m_ssest=ssest(data_id,nx);

figure
resid(m_ssest,data_v,'corr',5);
compare(m_ssest,data_g);

%%
H_ssest_th_dec=tf(m_ssest_th_dec);
H_ssest_th_dec=tf(3.935,[1 0]);
H_ssest_w_dec=tf(m_ssest_w_dec);
H_ssest=series(H_ssest_th_dec,H_ssest_w_dec);


