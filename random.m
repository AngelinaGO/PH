for k=1:100000
rand('state',1)
d5=rand;
plot(d5)
hold on
bb=atan(tan(Sun/1000)*sqrt(d5)); 
end