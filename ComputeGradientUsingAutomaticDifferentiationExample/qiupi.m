function [pi] = qiupi(k)
format long
l=1;
for n=1:k;
    OD=sqrt(1-(0.5*l)^2);
    CD=1-OD;
    l=sqrt(CD^2+(0.5*l)^2);
    pi=l*3*(2^n)
end
disp(pi)