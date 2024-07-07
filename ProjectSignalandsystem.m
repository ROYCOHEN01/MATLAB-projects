

%----- Question 1 ------ %
%{
t = linspace(0,1)
s = @(t) 1 - sqrt(t).*exp(-t)
plot(t,s(t))
t_1 =linspace(-1,0)
s_1 = @(t) -1 + sqrt(-t).*exp(t)
plot(t_1,s_1(t))


for i = -3:4
    
    if mod(i,2) == 0
        plot(linspace(i-1,i),s_1(t_1),'b')
    else
        plot(linspace(i-1,i),s(t),'b')

    end
    hold on
    
end
hold off

xlabel('t [second]')
ylabel('x(t)')

xlim([-4,4])
ylim([-1,1])




%this is sub question b

s_1 = @(t) -1 + sqrt(-t).*exp(t)
s = @(t) 1 - sqrt(t).*exp(-t)
n = 1:100
f_0_100 = -1 + (2*n-1)*2/400
p = pi
f_100_200 = -1 + (2*(n+100)-1)*2/400
k = 1:40

bk_p = s(f_100_200)*sin(p*(k .*(f_100_200.')))
bk_n = s_1(f_0_100)*sin(p*(k .*(f_0_100.')))
b_k = (bk_p+bk_n)/100

o = ones(40)
lower = tril(o);
x_m = @(t)  (lower .*b_k) * sin(p*((k) .* t.').');
%----------
s_mat = repmat(s(f_100_200),40,1)
e_x = ((x_m(f_100_200)-s_mat).^2)/200
s_mat1 = repmat(s_1(f_0_100),40,1)

e_x1 = (x_m(f_0_100)-s_mat1).^2/200
e_x_total = sum(e_x + e_x,2)


scatter(k,e_x_total)
xlabel('x_m')
ylabel('e_m')






x_matrix = x_m(linspace(-4,4,1000));

plot(linspace(-4,4,1000),x_matrix(30,:))

hold on

plot(linspace(-4,4,1000),x_matrix(15,:))

hold off

hold on

plot(linspace(-4,4,1000),x_matrix(5,:))

hold off
legend('x_3_0','x_1_5','x_5')

xlabel('t [seconds]')

ylabel('x_n(t)')
%}
%------- this is question 2 ----- %

x_c = @(t) exp(-((t-15).^2)/2)

samples = 1:51
x_n = x_c(samples)

%{
plot(linspace(0,80,1000),x_c(linspace(0,80,1000)))

hold on
stem(samples,x_n)
hold off

legend('x_c(t)','x[n]')
xlabel('t [second]')
%}

%------

%reiman_sum = sum(x_n.*exp(1:51))
%{
reiman_sum_F = @(x_n_w,TS) 0.3*sum(x_n_w.*exp(-1*TS))
T_n = 0.3*(1:281)


R_F = (real(Fourier)).^2
I_F = (imag(Fourier)).^2


plot(R_F)
hold on
plot ((R_F+I_F).^0.5)
hold off
%}

t_p =  5*(0:281)
x_n_c = x_c(t_p)

[w_s,Fourier] = FourierTransform(x_n_c,t_p)

display(Fourier)


R_F = (real(Fourier))
I_F = (imag(Fourier))

plot(w_s,R_F)
hold on
plot(w_s,I_F)
hold off

abs = (R_F.^2+I_F.^2).^0.5

phase = angle(Fourier)
display(w_s)
plot(w_s,abs)
hold off
hold on 
plot(w_s,phase)
hold off
xlabel('w')
phase_s = [char(8736),'F']
abs_s = '|F|'

legend(abs_s,phase_s)




%{
plot(t_p,x_n_c)
xlim([0 80])

hold on
plot(linspace(0,80,1000),x_c(linspace(0,80,1000)))
xlim([0 80])


hold off
%}

function [w_range,F] = FourierTransform(x_n,time_points)

    reiman = @(x_n_w,time_p, w_sum,TS) sum(TS*x_n_w .* exp(-time_p * i * w_sum));

    T_s = (time_points(:,2) - time_points(:,1));
    w = (2*pi)/T_s;
    w_range = linspace(-w/2,w/2,length(time_points));

    F_1 = zeros(1,length(w_range));
    
    for w_1 = 1:length(w_range)
        F_1(w_1) = reiman(x_n,time_points,w_range(w_1),T_s);
    end

    F = F_1;
    
end



