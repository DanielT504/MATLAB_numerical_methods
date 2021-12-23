% Q1

%setting up
clear;
clc;
Resistance=input("Q1\nEnter resistance (in Ohms): ");
iterations=[0 0];
ARAP=0;
syms T;

%assuming T between -200C and 850C, thus R between 18.52008 and 390.481125

%given resistance formulae:
if Resistance < 100
    f = 100*(1 + 3.9083*(10^(-3))*T - 5.775*(10^(-7))*T*T - 4.183*(10^(-12))*(T-100)*T*T*T) - Resistance;
    %worst case bounds
    L=-200;
    R=0;
else
    f = 100*(1 + 3.9083*(10^(-3))*T - 5.775*(10^(-7))*T*T) - Resistance;
    %worst case bounds
    L=0;
    R=850;
end


%bisection method:
x=0;
bound=[L R];
mid=(bound(1)+bound(2))/2;
oldMid=1000;

%iterate until appropriate error is reached
while x == 0
    x=1;
    iterations(1)=iterations(1)+1;

    mid=(bound(1)+bound(2))/2;
    
    ARAP=abs(100*(mid-oldMid)/mid);
    
    %determining if mid*left is above or below 0
    if subs(f,T,bound(1))*subs(f,T,mid) > 0
        bound(1)=mid;
        
        %end algorithm if appropriate error is reached
        if ARAP > 0.05
            x=0;
        end        
    elseif subs(f,T,bound(1))*subs(f,T,mid) < 0
        bound(2)=mid;
        
        %end algorithm if appropriate error is reached
        if ARAP > 0.05
            x=0;
        end        
    end
    
    oldMid=mid;
end


%Newton Raphson method:
Xi=(L+R)/2;
y=0;

%iterate until appropriate error is reached
while y == 0
    y=1;
    iterations(2)=iterations(2)+1;
    
    %formula for Newton's iteration
    Xj=Xi-double(subs(f,Xi)/subs(diff(f,T),Xi));
    
    ARAP(2)=abs(100*(Xj-Xi)/Xj);
    
    %exit when appropriate error is reached
    if ARAP(2) > 0.05
        y=0;
    end
    
    Xi=Xj;
end

%conclusion
fprintf("\nThe temperature obtained by bisection is %dC\n", mid);
fprintf("The temperature obtained by Newton Raphson is %dC\n", Xj);
fprintf("The number of required iterations for bisection is %i\n", iterations(1));
fprintf("The number of required iterations for Newton Raphson is %i\n", iterations(2));
fprintf("The absolute relative approximate error %% for bisection is %d%%\n", ARAP(1));
fprintf("The absolute relative approximate error %% for Newton Raphson is %d%%\n\n", ARAP(2));



% Q2

%reset variables and receive prompted input
clear;
syms x y z f1 f2 f3 f;
xi=input("\nQ2\nEnter initial value for x: ");
yi=input("Enter initial value for y: ");
zi=input("Enter initial value for z: ");
disp("*Note that coefficients must be expanded (ex. 2*y instead of 2y)");
f1=input("Enter equation 1:    0 = f1(x,y,z) = ");
f2=input("Enter equation 2:    0 = f2(x,y,z) = ");
f3=input("Enter equation 3:    0 = f3(x,y,z) = ");

%setting up
Xi=[xi yi zi];
Xj=[0 0 0];
f=[f1 f2 f3];
fColumn=[f1;f2;f3];
a=0;
iterations=0;

%I basically took my Newton Raphson method from above, and plugged in my
%Gauss elimination from Simulation 1 to find f(x,y,z)/f'(x,y,z)
while a == 0
    a=1;
    iterations=iterations+1;
    
    %creating an augmented matrix using the Jacobian as matrix A and
    %f1,f2,f3 as vector B, after substituting x,y,z values into both
    jac=jacobian(f,[x,y,z]);
    MyMatrixAug=[double(subs(jac,{x y z},{Xi(1) Xi(2) Xi(3)})) double(subs(fColumn,{x y z},{Xi(1) Xi(2) Xi(3)}))];
    
    %just Gauss elimination from Simulation 1:
    
    %forward elimination in each column
    for i=1:2
        afterI=MyMatrixAug([i:3],i);
        [val,row]=max(abs(afterI));

        row=row+i-1;

        %partial pivoting
        if row ~= i
            tmp1=MyMatrixAug(i,:);
            MyMatrixAug(i,:)=MyMatrixAug(row,:);
            MyMatrixAug(row,:)=tmp1;
        end

        %switching pivot row
        for j=(i+1):3
            tmp2=(MyMatrixAug(j,i)/MyMatrixAug(i,i))*MyMatrixAug(i,:);
            MyMatrixAug(j,:)=MyMatrixAug(j,:)-tmp2;
        end
    end

    Soln=MyMatrixAug(:,(4));

    %backward substitution for each row
    for m=1:3
        h=4-m;

        for k=(h+1):3
            Soln(h)=Soln(h)-MyMatrixAug(h,k)*Soln(k);
        end

        Soln(h)=Soln(h)/MyMatrixAug(h,h);
    end
    
    %formula for Newton's iteration
    Xj=Xi-Soln';
    
    ARAP=[abs(100*(Xj(1)-Xi(1))/Xj(1)) abs(100*(Xj(2)-Xi(2))/Xj(2)) abs(100*(Xj(3)-Xi(3))/Xj(3))];
    
    %exit when all three variables reach the appropriate error
    for h=1:3
        if ARAP(h) > 0.1
            a=0;
        end
    end
 
    Xi=Xj;
end

%conclusion
fprintf("\nThe values from your system of equations are as follows:\n");
fprintf("x = %d , y = %d , z = %d\n", Xi(1), Xi(2), Xi(3));
fprintf("The number of required iterations was %i\n", iterations);
fprintf("The absolute relative approximate error %% for x is %d%%\n", ARAP(1));
fprintf("The absolute relative approximate error %% for y is %d%%\n", ARAP(2));
fprintf("The absolute relative approximate error %% for z is %d%%\n", ARAP(3));