function [ K ] = solving_cracks( sigma,X1,X2 )

format longG

Length=X2-X1;

Center=(X1+X2)/2;

M=zeros(length(Center));

I=eye(length(Center));


for i=1:length(Center)  %sending crack
    
   for j=1:length(Center)  %receiving crack
    
         if abs(i-j)>0
             

            if i<j
                djk=X2(j)-Center(i);
                cjk=X1(j)-Center(i);
                ak=Length(i)/2;
            else
                djk=Center(i)-X1(j);
                cjk=Center(i)-X2(j);
                ak=Length(i)/2;   
            end
                
            M(j,i)=1/(djk-cjk) * (cjk-djk-(cjk^2-ak^2)^0.5+(djk^2-ak^2)^0.5);
         
         end

   end
   
end
'M';
M;
(I-M);

'p bar';
P_bracket=(I-M)\ones(length(Center),1)*(sigma);    %average traction for each crack


m=10;     %integration points in a crack
Integral=[0.2093662482094422,0.9991601288170098
0.4642615919622012,0.9776208824732407
0.6455583892806666,0.8758595017700398
0.6891745658564888,0.6293961480326324
0.5717773617691598,0.231725670698451
0.3566393715525271,-0.231725670698451
0.1567517813840062,-0.629396148032632
4.272171770125290E-02,-0.87585950177004
5.253668598345737E-03,-0.977620882473241
8.795727567153020E-05,-0.99916012881701];

Points=Integral(:,2)/2+0.5;
Weights=Integral(:,1);

P=zeros(length(Center),m);
B=zeros(length(Center),m);

for j=1:length(Center)
    for i=1:m
        B(j,i)=X1(j)*(1-Points(i))+X2(j)*Points(i);
 
    end
end

B;

for k=1:length(Center)  % receiver
    for j=1:length(Center)  % sender
        
        if abs(j-k)>0
            for i=1:m
                
                
                r=abs(B(k,i)-Center(j));
                r1=abs(B(k,i)-X1(j));
                r2=abs(B(k,i)-X2(j));
                P(k,i)=P(k,i)+(r/(r1*r2)^0.5-1)   *P_bracket(j);
                
            end
        end
        
        
    end
    
    for i=1:m
        P(k,i)=P(k,i)+sigma;

    end
end

%  P(k,i) is the stress on k crack at the ith points
P;
sum(P,2);

K=[];

for j=1:length(Center)
    K1=0;
    K2=0;
    a=Length(j)/2;
    
    for i=1:m

        K1=K1+P(j,i)*Weights(m-i+1);
        K2=K2+P(j,i)*Weights(i);
        
    end
    K(j,1)=K1*(a/pi)^0.5;
    K(j,2)=K2*(a/pi)^0.5;

end
K;


end

