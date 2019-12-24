%% 定义
Sun=4.65;
N_W=10000;
N_L=1;
N_Num=100;
M=90;
D1=0.064;
D2=0.070;
D3=0.116;
D4=0.120;
R2=D2/2;
R4=D4/2;
DNI=1367;
pi=3.1415926;

L=150;
f=1.71;
W=5.76;
reflectivity=0.95;
transmittance=0.95;
Glass_reflect=0.03;
absorption=1.00;
X=W/2;
Dx=W/N_W;
Da=180/M;

dd=0;   %入射光线与垂直方向的夹角
%% 光线数量统计
    for i=1:M
    Num_R(i,1)=0; 
    L1_Num(i,1)=0;
    R1_Num(i,1)=0;
    Num_L(i,1)=0; 
    R2_Num(i,1)=0;
    L2_Num(i,1)=0;
    Flux_R(i,1)=0; 
    Flux_L(i,1)=0; 
    end

 h=waitbar(0);
%% 产生随机数
for k=1:N_W
    rand('state',1)
       d1=rand;
       d2=rand;
       d3=rand;
       d4=rand;
       d5=rand;
       d6=rand;
       d9=rand;
       d10=rand;
       d11=rand;
       
       x0=Dx/2+(k-N_W/2)*Dx; 
	   y0=x0*x0/(4*f); 
       plot(x0,y0,'b')
       hold on
       
       
       if x0>-R2&&x0<0               %直射至金属管左侧          
           thetal=-180-asin(x0/R2)*180/pi;
           for n=1:M
               s1=-(n-1)*Da;
               s2=-n*Da;
               if  thetal>=s2&&thetal<=s1
                    L2_Num(n)=L2_Num(n)+N_Num*transmittance*absorption;
                   break
               end
           end
       elseif  x0>0&&x0<R2           %直射至金属管右侧
           thetal=180-asin(x0/R2)*180/pi;
           for n=1:M
               s1=(n-1)*Da;
               s2=n*Da;
                if  thetal>=s1&&thetal<=s2
                    R1_Num(n)=R1_Num(n)+N_Num*transmittance*absorption;
                    break
                end
           end     
       elseif x0<=-R2&&x0>-R4         %直射至左侧真空范围
           for i=1:N_Num
               if d3<=transmittance
                   aa=2*pi*d4;         %光锥的方向角
                  bb=atan(tan(Sun/1000)*sqrt(d5));     %光锥的顶角
                   
                   x111=sin(bb)*cos(aa);     %光锥内的光线，以底面圆为坐标轴
                   y111=cos(bb);
                   z111=sin(bb)*sin(aa);
                   
                   x1=x111;
                   y1=cos(dd)*y111+sin(dd)*z111;
                   z1=-sin(dd)*y111+cos(dd)*z111;
                   
                   if d6<=Glass_reflect    
                       Nx0=-x0;
                       Ny0=sqrt(R4*R4-x0*x0);
                       Nz0=0;
                       y0=f-Ny0;
                       P0N0=sqrt(Nx0*Nx0+Ny0*Ny0+Nz0*Nz0);
                       S=2*(x1*Nx0+y1*Ny0+z1*Nz0)/(P0N0*P0N0);
                       x2=S*Nx0-x1;
                       y2=S*Ny0-y1;
                       z2=S*Nz0-z1;
                       
                       E=x2*x2+y2*y2;
			           F=2*x2*y2*(y0-f)-2*y2*y2*x0;
			           G=y2*y2*x0*x0+x2*x2*(y0-f)*(y0-f)-2*x2*y2*x0*(y0-f)-x2*x2*R2*R2;
			           D=F*F-4*E*G;
                       if D>=0
                           x3=(-F-sqrt(D))/(2*E);
                           y3=y2/x2*(x3-x0)+y0;
                           if d9<=absorption      %光线经玻璃管内侧反射被吸热管吸收
                               if x3<0
                                   if y3<=f
                                       thetal=asin(x3/R2)*180/pi;
                                   elseif y3>f
                                       thetal=(-asin(x3/R2)-pi)*180/pi;
                                   end
                               elseif x3>=0
                                   if y3<=f
                                       thetal=asin(x3/R2)*180/pi;
                                   elseif y3>f
                                       thetal=(pi-asin(x3/R2))*180/pi;
                                   end
                               end
                           end
                           for n=1:M
                               s1=-(n-1)*Da;
				               s2=-n*Da;
				               if thetal>=s2&&thetal<s1
				                   L2_Num(n)=L2_Num(n)+1; 
					                break
                               end
                           end
                           for n=1:M
                               s1=(n-1)*Da;
				               s2=n*Da;
				               if thetal>=s1&&thetal<s2
				                      R2_Num(n)=R2_Num(n)+1;
					                break
                               end
                           end
                       end
                   elseif d6>Glass_reflect&&d6<=Glass_reflect+transmittance        %光线到达抛物面反射镜
                      Nx0=-x0;
		              Ny0=2*f;
		              Nz0=0;
                      x1_0=0;                 
		              y1_0=1;
		              z1_0=0;
		              P0N0=sqrt(Nx0*Nx0+Ny0*Ny0+Nz0*Nz0);
		              P0N1=sqrt(x1_0*x1_0+y1_0*y1_0+z1_0*z1_0);
                      cc =2*acos((Nx0*x1_0+Ny0*y1_0+Nz0*z1_0)/(P0N0*P0N1));
                      if(d11<=reflectivity)
                          x2=-cos(-cc)*sin(bb)*cos(aa)-sin(-cc)*cos(dd)*cos(bb)-sin(-cc)*sin(dd)*sin(bb)*sin(aa); 
                          y2=-sin(-cc)*sin(bb)*cos(aa)+cos(-cc)*cos(dd)*cos(bb)+cos(-cc)*sin(dd)*sin(bb)*sin(aa);
                          z2=sin(dd)*cos(bb)-cos(dd)*sin(bb)*sin(aa);
                           
                          E=x2*x2+y2*y2;
			              F=2*x2*y2*(y0-f)-2*y2*y2*x0;
			              G=y2*y2*x0*x0+x2*x2*(y0-f)*(y0-f)-2*x2*y2*x0*(y0-f)-x2*x2*R2*R2;
			              D=F*F-4*E*G;
                          if D>=0
                           x3=(-F-sqrt(D))/(2*E);
                           y3=y2/x2*(x3-x0)+y0;
                           if d10<=transmittance
                             if d9<=absorption      %光线经玻璃管内侧反射被吸热管吸收
                               if x3<0
                                   if y3<=f
                                       thetal=asin(x3/R2)*180/pi;
                                   elseif y3>f
                                       thetal=(-asin(x3/R2)-pi)*180/pi;
                                   end
                               elseif x3>=0
                                   if y3<=f
                                       thetal=asin(x3/R2)*180/pi;
                                   elseif y3>f
                                       thetal=(pi-asin(x3/R2))*180/pi;
                                   end
                               end
                             end
                           for n=1:M
                               s1=-(n-1)*Da;
				               s2=-n*Da;
				               if thetal>=s2&&thetal<s1
				                      L2_Num(n)=L2_Num(n)+1; 
					                break
                                end
                           end
                           for n=1:M
                               s1=(n-1)*Da;
				               s2=n*Da;
				               if thetal>=s1&&thetal<s2
				                      R2_Num(n)=R2_Num(n)+1;
					                break
                               end
                           end
                           end
                          end
                      end
                 end   
             end
           end
           
       elseif x0>=R2&&x0<R4      %直射至右侧真空范围
           for i=1:N_Num
               if d3<=transmittance
                   aa=2*pi*d4;         %光锥的方向角
                   bb=atan(tan(Sun/1000)*sqrt(d5));     %光锥的顶角
                   
                   x111=sin(bb)*cos(aa);     %光锥内的光线，以底面圆为坐标轴
                   y111=cos(bb);
                   z111=sin(bb)*sin(aa);
                   
                   x1=x111;
                   y1=cos(dd)*y111+sin(dd)*z111;
                   z1=-sin(dd)*y111+cos(dd)*z111;
                   
                   if d6<=Glass_reflect    
                       Nx0=-x0;
                       Ny0=sqrt(R4*R4-x0*x0);
                       Nz0=0;
                       y0=f-Ny0;
                       P0N0=sqrt(Nx0*Nx0+Ny0*Ny0+Nz0*Nz0);
                       S=2*(x1*Nx0+y1*Ny0+z1*Nz0)/(P0N0*P0N0);
                       x2=S*Nx0-x1;
                       y2=S*Ny0-y1;
                       z2=S*Nz0-z1;
                       
                       E=x2*x2+y2*y2;
			           F=2*x2*y2*(y0-f)-2*y2*y2*x0;
			           G=y2*y2*x0*x0+x2*x2*(y0-f)*(y0-f)-2*x2*y2*x0*(y0-f)-x2*x2*R2*R2;
			           D=F*F-4*E*G;
                       if D>=0
                           x3=(-F-sqrt(D))/(2*E);
                           y3=y2/x2*(x3-x0)+y0;
                           if d9<=absorption      %光线经玻璃管内侧反射被吸热管吸收
                               if x3<0
                                   if y3<=f
                                       thetal=asin(x3/R2)*180/pi;
                                   elseif y3>f
                                       thetal=(-asin(x3/R2)-pi)*180/pi;
                                   end
                               elseif x3>=0
                                   if y3<=f
                                       thetal=asin(x3/R2)*180/pi;
                                   elseif y3>f
                                       thetal=(pi-asin(x3/R2))*180/pi;
                                   end
                               end
                           end
                           for n=1:M
                               s1=-(n-1)*Da;
				               s2=-n*Da;
				               if thetal>=s2&&thetal<s1
				                   L1_Num(n)=L1_Num(n)+1; 
					                break
                               end
                           end
                           for n=1:M
                               s1=(n-1)*Da;
				               s2=n*Da;
				               if thetal>=s1&&thetal<s2
				                      R1_Num(n)=R1_Num(n)+1;
					                break
                               end
                           end
                       end
                   elseif d6>Glass_reflect&&d6<=Glass_reflect+transmittance        %光线到达抛物面反射镜
                      Nx0=-x0;
		              Ny0=2*f;
		              Nz0=0;
                      x1_0=0;                 
		              y1_0=1;
		              z1_0=0;
		              P0N0=sqrt(Nx0*Nx0+Ny0*Ny0+Nz0*Nz0);
		              P0N1=sqrt(x1_0*x1_0+y1_0*y1_0+z1_0*z1_0);
                      cc =2*acos((Nx0*x1_0+Ny0*y1_0+Nz0*z1_0)/(P0N0*P0N1));
                      if d11<=reflectivity
                          x2=cos(cc)*sin(bb)*cos(aa)-sin(cc)*cos(dd)*cos(bb)-sin(cc)*sin(dd)*sin(bb)*sin(aa); 
                          y2=-sin(cc)*sin(bb)*cos(aa)+cos(cc)*cos(dd)*cos(bb)+cos(cc)*sin(dd)*sin(bb)*sin(aa);
                          z2=sin(dd)*cos(bb)-cos(dd)*sin(bb)*sin(aa);
                           
                          E=x2*x2+y2*y2;
			              F=2*x2*y2*(y0-f)-2*y2*y2*x0;
			              G=y2*y2*x0*x0+x2*x2*(y0-f)*(y0-f)-2*x2*y2*x0*(y0-f)-x2*x2*R2*R2;
			              D=F*F-4*E*G;
                          if D>=0
                           x3=(-F-sqrt(D))/(2*E);
                           y3=y2/x2*(x3-x0)+y0;
                           if d10<=transmittance
                             if d9<=absorption      %光线经玻璃管内侧反射被吸热管吸收
                               if x3<0
                                   if y3<=f
                                       thetal=asin(x3/R2)*180/pi;
                                   elseif y3>f
                                       thetal=(-asin(x3/R2)-pi)*180/pi;
                                   end
                               elseif x3>=0
                                   if y3<=f
                                       thetal=asin(x3/R2)*180/pi;
                                   elseif y3>f
                                       thetal=(pi-asin(x3/R2))*180/pi;
                                   end
                               end
                             end
                           for n=1:M
                               s1=-(n-1)*Da;
				               s2=-n*Da;
				               if thetal>=s2&&thetal<s1
				                      L1_Num(n)=L1_Num(n)+1; 
					                break
                                end
                           end
                           for n=1:M
                               s1=(n-1)*Da;
				               s2=n*Da;
				               if thetal>=s1&&thetal<s2
				                      R1_Num(n)=R1_Num(n)+1;
					                break
                               end
                           end
                           end
                          end
                      end
                 end   
             end
           end
               
       elseif x0<=-R4&&x0>=-X          %直射在左侧反射镜上
           Nx0=-x0;
		   Ny0=2*f;
		   Nz0=0;
           x1_0=0;                 
		   y1_0=1;
		   z1_0=0;
           P0N0=sqrt(Nx0*Nx0+Ny0*Ny0+Nz0*Nz0);
		   P0N1=sqrt(x1_0*x1_0+y1_0*y1_0+z1_0*z1_0);
           
           cc=2*acos((Nx0*x1_0+Ny0*y1_0+Nz0*z1_0)/(P0N0*P0N1));
           for i=1:N_Num
               if d11<=reflectivity
                   aa=2*pi*d1;         %光锥的方向角
                   bb=atan(tan(Sun/1000)*sqrt(d5));
                   
                   x2=-cos(-cc)*sin(bb)*cos(aa)-sin(-cc)*cos(dd)*cos(bb)-sin(-cc)*sin(dd)*sin(bb)*sin(aa);
			       y2=-sin(-cc)*sin(bb)*cos(aa)+cos(-cc)*cos(dd)*cos(bb)+cos(-cc)*sin(dd)*sin(bb)*sin(aa);
			       z2=+sin(dd)*cos(bb)-cos(dd)*sin(bb)*sin(aa);
                   E=x2*x2+y2*y2;
			       F=2*x2*y2*(y0-f)-2*y2*y2*x0;
			       G=y2*y2*x0*x0+x2*x2*(y0-f)*(y0-f)-2*x2*y2*x0*(y0-f)-x2*x2*R2*R2;
			       D=F*F-4*E*G;
                   if D>=0
                           x3=(-F-sqrt(D))/(2*E);
                           y3=y2/x2*(x3-x0)+y0;
                           if d10<=transmittance
                             if d9<=absorption      %光线经玻璃管内侧反射被吸热管吸收
                               if x3<0
                                   bb=bb*1000;
                                   if y3<=f
                                       thetal=asin(x3/R2)*180/pi;
                                   elseif y3>f
                                       thetal=(-asin(x3/R2)-pi)*180/pi;
                                   end
                               elseif x3>=0
                                   if y3<=f
                                       thetal=asin(x3/R2)*180/pi;
                                   elseif y3>f
                                       thetal=(pi-asin(x3/R2))*180/pi;
                                   end
                               end
                             end
                           for n=1:M
                               s1=-(n-1)*Da;
				               s2=-n*Da;
				               if thetal>=s2&&thetal<s1
				                      L2_Num(n)=L2_Num(n)+1; 
					                break
                                end
                           end
                           for n=1:M
                               s1=(n-1)*Da;
				               s2=n*Da;
				               if thetal>=s1&&thetal<s2
				                      R2_Num(n)=R2_Num(n)+1;
					                break
                               end
                           end
                           end
                   end
               end
           end
       elseif x0>=R4&&x0<=X
            Nx0=-x0;
			Ny0=2*f;
			Nz0=0;
            
             x1_0=0;                 
			 y1_0=1;
			 z1_0=0;
             P0N0=sqrt(Nx0*Nx0+Ny0*Ny0+Nz0*Nz0);
			 P0N1=sqrt(x1_0*x1_0+y1_0*y1_0+z1_0*z1_0);
             cc=2*acos((Nx0*x1_0+Ny0*y1_0+Nz0*z1_0)/(P0N0*P0N1));
             for i=1:N_Num
               if d11<=reflectivity
                   aa=2*pi*d1;         %光锥的方向角
                   bb=atan(tan(Sun/1000)*sqrt(d2));
                   x2=-cos(cc)*sin(bb)*cos(aa)-sin(cc)*cos(dd)*cos(bb)-sin(-cc)*sin(dd)*sin(bb)*sin(aa);
			       y2=-sin(cc)*sin(bb)*cos(aa)+cos(cc)*cos(dd)*cos(bb)+cos(cc)*sin(dd)*sin(bb)*sin(aa);
			       z2=+sin(dd)*cos(bb)-cos(dd)*sin(bb)*sin(aa);
                   E=x2*x2+y2*y2;
			       F=2*x2*y2*(y0-f)-2*y2*y2*x0;
			       G=y2*y2*x0*x0+x2*x2*(y0-f)*(y0-f)-2*x2*y2*x0*(y0-f)-x2*x2*R2*R2;
			       D=F*F-4*E*G;
                   if D>=0
                           x3=(-F-sqrt(D))/(2*E);
                           y3=y2/x2*(x3-x0)+y0;
                           if d10<=transmittance
                             if d9<=absorption      %光线经玻璃管内侧反射被吸热管吸收
                               if x3<0
                                   bb=bb*1000;
                                   if y3<=f
                                       thetal=asin(x3/R2)*180/pi;
                                   elseif y3>f
                                       thetal=(-asin(x3/R2)-pi)*180/pi;
                                   end
                               elseif x3>=0
                                   if y3<=f
                                       thetal=asin(x3/R2)*180/pi;
                                   elseif y3>f
                                       thetal=(pi-asin(x3/R2))*180/pi;
                                   end
                               end
                             end
                           for n=1:M
                               s1=-(n-1)*Da;
				               s2=-n*Da;
				               if thetal>=s2&&thetal<s1
				                      L1_Num(n)=L1_Num(n)+1; 
					                break
                                end
                           end
                           for n=1:M
                               s1=(n-1)*Da;
				               s2=n*Da;
				               if thetal>=s1&&thetal<s2
				                      R1_Num(n)=R1_Num(n)+1;
					                break
                               end
                           end
                           end
                   end
               end
             end
       end

 str=[num2str(k/N_W),'2%'];
 waitbar(k/N_W,h,str)
 %waitbar(k/N_W);
 
      
       
end
       for j=1:M
          Num_R(j)=R1_Num(j)+R2_Num(j); 
       end
       for j=1:M
          Num_L(j)=L1_Num(j)+L2_Num(j); 
       end
  
   alpha=0:pi/20:2*pi;
   X_abs=0.035*cos(alpha);
   Y_abs=0.035*sin(alpha)+1.71;
   plot(X_abs,Y_abs,'c')
   hold on

    Q=DNI*W;
    q1=Q/(N_W*N_L*N_Num);
    ss=pi*R2/M;
    Ave_Flux=0;
   
    for i=1:M
	  Flux_R(i)=q1*Num_R(i)/ss;
	  Gr_R=Flux_R(i)/DNI;
      sprintf('%.2f',Gr_R)
	  Ave_Flux=Ave_Flux+Flux_R(i);
    end
    for i=1:M
	  Flux_L(i)=q1*Num_L(i)/ss;
	  Gr_L=Flux_L(i)/DNI;
      sprintf('%.2f',Gr_L)
	  Ave_Flux=Ave_Flux+Flux_L(i);
    
    end
    
    
    Ave_Flux=Ave_Flux/(2*M);
    
    Total_Energy=Ave_Flux*(2*pi*R2)/transmittance;
    Glass_Absor=Total_Energy*(1-transmittance-Glass_reflect);
    Glass_Flux=Glass_Absor/(2*pi*R4);

    sprintf('吸热管平均能流密度%.8f',Ave_Flux)
    sprintf('玻璃管吸收的平均能流密度%.8f',Glass_Flux)
    
    


        
 
 
 