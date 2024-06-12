clear
clc

format longG

L=1; 

Stress=1;  
Kc=Stress*(pi*L/2)^0.5;

num_L=100;
sample_size=L*num_L;

number_of_samples=20000;


N=num_L+1;  %num of nodes


AVE_S=[];
STD_S=[];
AVE_S_ini=[];
STD_S_ini=[];


for designed_ratio=0.1:0.1:0.9   
    
    level=norminv(designed_ratio,0,0.81);


for ith_sample=0:number_of_samples-1
    
    ith_sample


Void=[];
void_count=0;

surf=normrnd(0,1,N,1);
surf=surf-level;
high_or_low=surf./(abs(surf));

void_start=1;
void_end=1;
void=[];

for i=2:length(surf)

    if high_or_low(i-1)<0 && high_or_low(i)>0 
        
        void_end=i-abs(surf(i))/(abs(surf(i))+abs(surf(i-1)));
        void=[void_start;void_end];
        void_count=void_count+1;
        Void(:,void_count)=void;
     
    elseif high_or_low(i-1)>0 && high_or_low(i)<0
        
        void_start=i-abs(surf(i))/(abs(surf(i))+abs(surf(i-1)));  
            
    elseif high_or_low(i)<0 && i==length(surf)
        
        void_end=length(surf);
        void=[void_start;void_end];
        void_count=void_count+1;
        Void(:,void_count)=void;
        
    end
    
end
    

    Void_X1=Void(1,:)*L;
    Void_X2=Void(2,:)*L;
    
    
    s=[];
    for j=1:1                 %times of crack expansion
        j;
        C_store=Void_X2-Void_X1;
        
        if isempty(C_store)==1
            break
        end
        
        K=solving_cracks( Stress,Void_X1,Void_X2 );
        K1=K(:,1);
        K2=K(:,2);
        
        K1_max=max(K1);
        K2_max=max(K2);
        K_max=max(K1_max,K2_max);
        
        if length(Void_X1)>=3
            if K1_max>K2_max && find(K1==K1_max)-1>0
                Void_X1(find(K1==K1_max))=[];
                Void_X2(find(K1==K1_max)-1)=[];
            elseif K1_max<=K2_max && find(K2==K2_max)+1<=length(Void_X1)
                Void_X2(find(K2==K2_max))=[];
                Void_X1(find(K2==K2_max)+1)=[];
            end
            
        else
            s(j)=Kc/K_max;
            break
        end
        
        %%%%
        
        s(j)=Kc/K_max;
        
        length_critical= (Kc/max(s))^2/pi*2;   
        Length_store=Void_X2-Void_X1;
        
        if max(Length_store) > length_critical
            j;
            break
        end
        
    end
      
    if isempty(s)==1
        s=9999
    end
    
    if isreal(s(1))==0 || isreal(max(s))==0
        return
    end
    
    S_initial(ith_sample+1)=s(1);        % initial stress
    S(ith_sample+1)=max(s);              % max stress
    
end


for i=1:20
    ave_S=mean(S);
    std_S=std(S);
    j=1;
    while j<=length(S)
        if S(j)>ave_S+6*std_S || S(j)<ave_S-6*std_S || S(j)==9999
            S(j)=[];
            S_initial(j)=[];
            j=j-1;
        end
        j=j+1;
    end
end

AVE_S=[AVE_S;mean(S)];
STD_S=[STD_S;std(S)];
AVE_S_ini=[AVE_S_ini;mean(S_initial)];
STD_S_ini=[STD_S_ini;std(S_initial)];


end

AVE_S
STD_S
AVE_S_ini
STD_S_ini


