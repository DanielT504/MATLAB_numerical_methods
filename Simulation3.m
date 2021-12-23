%setting up
clc;
clear;

%loading data
MyMatrix=load('test1.txt');
%MyMatrix=load('test2.txt');
Size=size(MyMatrix,1);
%intro
disp("Select the function to fit your data (1, 2, or 3):");
disp("1. Polynomial: (y=(a0)+(a1)x+...+(am)x^m)");
disp("2. Exponential: (y=ae^(bx))");
disp("3. Saturation: (y=(ax)/(b+x))");
which=input("Choice: ");

x=MyMatrix(:,1);
y=MyMatrix(:,2);

%warning cases
zeroDiv=0;
zeroVal=0;
neg=0;
for i=1:Size
    if y(i) < 0
        neg=1;
    end
    
    if y(i) == 0 || x(i) == 0
        zeroVal=1;
    end
end

%render increments of 0.01 in figure
x0=min(MyMatrix(:,1)):0.01:max(MyMatrix(:,1));

%polynomial regression
if which == 1
    deg=input("What degree (1, 2, or 3): ");
    
    a0=0;
    a1=0;
    a2=0;
    a3=0;
    
    %set to x just to have the appropriate size
    second=x;
    third=x;
    fourth=x;
    fifth=x;
    sixth=x;
    product=x;
    productx2=x;
    productx3=x;

    %each series (later to be summed) required in the regression matrix
    for i=1:Size
        second(i)=(x(i))^2;
        third(i)=(x(i))^3;
        fourth(i)=(x(i))^4;
        fifth(i)=(x(i))^5;
        sixth(i)=(x(i))^6;
        product(i)=y(i)*x(i);
        productx2(i)=y(i)*((x(i))^2);
        productx3(i)=y(i)*((x(i))^3);
    end
    
    %linear
    if deg == 1
        %check for division by zero
        if (Size*sum(second)-(sum(x))^2) ~= 0
            %linear constant formula
            a1=(Size*sum(product)-sum(x)*sum(y))/(Size*sum(second)-(sum(x))^2);
        else
            zeroDiv = 1;
        end
        
        %linear constant formula
        a0=mean(y)-a1*mean(x);
    %quadratic
    elseif deg == 2
        %setting regression matrix and vector
        Matrix=[Size sum(x) sum(second); sum(x) sum(second) sum(third); sum(second) sum(third) sum(fourth)];
        Vector=[sum(y); sum(product); sum(productx2)];
        a0a1a2=(inv(Matrix))*Vector; %simple matrix inversion
        
        %determining constants
        a0=a0a1a2(1);
        a1=a0a1a2(2);
        a2=a0a1a2(3);
    %cubic
    else
        Matrix=[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
        Matrix(1,:)=[Size sum(x) sum(second) sum(third)];            %
        Matrix(2,:)=[sum(x) sum(second) sum(third) sum(fourth)];     % could have been done in one line
        Matrix(3,:)=[sum(second) sum(third) sum(fourth) sum(fifth)]; % (just setting regression matrix)
        Matrix(4,:)=[sum(third) sum(fourth) sum(fifth) sum(sixth)];  %
        Vector=[sum(y); sum(product); sum(productx2); sum(productx3)]; % regression vector
        a0a1a2a3=(inv(Matrix))*Vector;  %simple matrix inversion
        
        %determining constants
        a0=a0a1a2a3(1);
        a1=a0a1a2a3(2);
        a2=a0a1a2a3(3);
        a3=a0a1a2a3(4);
    end

    if zeroDiv == 0
        % St and Sr
        StVec=x;
        SrVec=x;

        %given formulae
        for i=1:Size
            StVec(i)=(y(i)-mean(y))^2;
            SrVec(i)=(y(i)-a0-a1*x(i)-a2*(x(i))^2-a3*(x(i))^3)^2;
        end
        
        %unless a division by zero will occur
        if sum(StVec) ~= 0
            %R^2 value
            R1=(sum(StVec)-sum(SrVec))/sum(StVec);

            %equation and title to plot
            Y=a0+x0.*a1+(x0.^2).*a2+(x0.^3).*a3;
            Title=['Polynomial, y=' num2str(a0) '+' num2str(a1) 'x+' num2str(a2) 'x^{2}+' num2str(a3) 'x^{3}, R^{2}=' num2str(R1)];
        else
            zeroDiv=1;
        end
    end
%exponential regression
elseif which == 2
    %check if a negative value will cause an invalid logarithm
    if neg == 0
        %each series (later to be summed) required for regression formulae
        linY=log(y);
        product=x;
        square=x;
        squareLnY=x;

        for i=1:Size
            product(i)=x(i)*linY(i);
            square(i)=(x(i))^2;
            squareLnY(i)=(linY(i))^2;
        end

        %simplified formulae
        Sxy=sum(product)-Size*mean(x)*mean(linY);
        Sxx=sum(square)-Size*(mean(x))^2;
        Syy=sum(squareLnY)-Size*(mean(linY))^2;
        
        %check for division by zero
        if Sxx ~= 0 && Syy ~= 0
            %determine constants
            b=Sxy/Sxx;
            a=exp(mean(linY)-b*mean(x));

            %R^2 value
            R2=(Sxy^2)/(Sxx*Syy);

            %equation and title to plot
            Y=a.*exp(x0.*b);
            Title=['Exponential, y=' num2str(a) 'e^{' num2str(b) 'x}, R^{2}=' num2str(R2)];
        else
            zeroDiv=1;
        end
    end
%saturation-growth-rate regression
else
    %check if a zero value will cause a division by zero
    if zeroVal == 0
        %each series (later to be summed) required for regression formulae
        invY=x;
        invX=x;
        product=x;
        square=x;

        for i=1:Size
            invY(i)=1/y(i);
            invX(i)=1/x(i);
            product(i)=invX(i)*invY(i);
            square(i)=(invX(i))^2;
        end

        %check for division by zero
        if (Size*sum(square)-(sum(invX))^2) ~= 0
            %given formulae for constants
            A1=(Size*sum(product)-sum(invX)*sum(invY))/(Size*sum(square)-(sum(invX))^2);
            A0=mean(invY)-A1*mean(invX);
            %linearizing constants
            A3=1/A0;
            B3=A1*A3;
            
            % St and Sr
            StVec=x;
            SrVec=x;

            %given formulae
            for i=1:Size
                StVec(i)=(y(i)-mean(y))^2;
                SrVec(i)=(y(i)-(A3*x(i))/(B3+x(i)))^2;
            end

            %unless a division by zero will occur
            if sum(StVec) ~= 0
                %R^2 value
                R3=(sum(StVec)-sum(SrVec))/sum(StVec);

                %equation and title to plot
                Y=A3.*x0./(B3+x0);
                Title=['Saturation, y=(' num2str(A3) 'x)/(' num2str(B3) '+x), R^{2}=' num2str(R3)];
            else
                zeroDiv=1;
            end
        else
            zeroDiv=1;
        end
    else
        zeroDiv=1;
    end
end

%plot if no warnings occurred
if (neg == 0 || which ~= 2) && zeroDiv == 0
    figure;
    plot(x,y,'.'); %raw data
    hold on;
    plot(x0,Y); %estimated regression line
    %axis to fit min and max data
    axis([min(MyMatrix(:,1)) max(MyMatrix(:,1)) min(MyMatrix(:,2)) max(MyMatrix(:,2))]);
    xlabel('x');
    ylabel('y');
    %predetermined title
    title(Title);
else
    %warnings
    if zeroDiv ~= 0
        disp("Warning: Division by zero");
    else
        disp("Warning: Invalid natural logarithm");
    end
end