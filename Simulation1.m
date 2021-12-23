clear;
clc;

%setting up
MyMatrix=load('A.txt'); %should be nxn
MyVector=load('B.txt'); %should be n values
Size=size(MyVector);
n=Size(:,1);



% Part a)

MyMatrixInv=inv(MyMatrix);
SolnA=MyMatrixInv*MyVector;
SolnARounded=round(SolnA,5,'significant');
disp("Mesh currents in increasing order (calculated using matrix inversion to 5 sig figs):");
disp(SolnARounded);



% Part b)

MyMatrixAug=[MyMatrix MyVector]; %concatenating B onto A to make an augmented matrix

%forward elimination in each column
for i=1:(n-1)
    afterI=MyMatrixAug([i:n],i);
    [val,row]=max(abs(afterI));
    
    row=row+i-1;
    
    %partial pivoting
    if row ~= i
        tmp1=MyMatrixAug(i,:);
        MyMatrixAug(i,:)=MyMatrixAug(row,:);
        MyMatrixAug(row,:)=tmp1;
    end
    
    %switching pivot row
    for j=(i+1):n
        tmp2=(MyMatrixAug(j,i)/MyMatrixAug(i,i))*MyMatrixAug(i,:);
        MyMatrixAug(j,:)=MyMatrixAug(j,:)-tmp2;
    end
end

SolnB=MyMatrixAug(:,(n+1));

%backward substitution for each row
for m=1:n
    h=(n+1)-m;
    
    for k=(h+1):n
        SolnB(h)=SolnB(h)-MyMatrixAug(h,k)*SolnB(k);
    end
    
    SolnB(h)=SolnB(h)/MyMatrixAug(h,h);
end

SolnBRounded=round(SolnB,5,'significant');
disp("Mesh currents in increasing order (calculated using Gaussian elimination to 5 sig figs):");
disp(SolnBRounded);



% Part c)

Error=[1 0.5 0.1 0.01];
%all four levels of accuracy
for s=1:4
    Iterations=0;
    SolnC=MyVector*0; %initial guess is all zeros

    ARAP=MyVector;
    for q=1:n
        ARAP(q)=100;
    end

    x=0;

    while x == 0
        SolnCTemp=SolnC;
        x=1;
        
        %rearranging equations
        for o=1:n
            sum=0;

            for p=1:n
                sum=sum+(SolnC(p)*MyMatrix(o,p));
            end
            sum=sum-SolnC(o)*MyMatrix(o,o);

            SolnC(o)=(MyVector(o)-sum)/MyMatrix(o,o);
        end

        %checking if all values have reached below the appropriate error
        for r=1:n
            ARAP(r)=100*abs((SolnC(r)-SolnCTemp(r))/SolnC(r));

            if ARAP(r) >= Error(s)
                x=0;
            end
        end
        
        Iterations=Iterations+1;
    end
    
    %output with requested information
    SolnCRounded=round(SolnC,5,'significant');
    disp("Mesh currents in increasing order (calculated using the Gauss-Seidel method to 5 sig figs)");
    disp("Stopped when absolute relative aproximate error reached (in %):");
    disp(Error(s));
    disp(SolnCRounded);
    disp("Iterations  needed to reach this accuracy:");
    disp(Iterations);
end
