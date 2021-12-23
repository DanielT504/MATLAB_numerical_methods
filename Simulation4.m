%setting up
clc;
clear;

%loading data
MyMatrix=load('test_1.txt');
Size=size(MyMatrix,1);
%intro
which=input("Derivative (1) or integral (2) : ");

x=MyMatrix(:,1);
y=MyMatrix(:,2);

%render increments of 0.01 in figure
x0=min(MyMatrix(:,1)):0.01:max(MyMatrix(:,1));

eqSpacing = 1;
minSeg = x(2)-x(1);

%checking for equal spacing and finding minimum space
if Size > 2
    for i=3:Size
        if x(i)-x(i-1) ~= x(2)-x(1)
            if x(i)-x(i-1) < minSeg
                minSeg = x(i)-x(i-1);
            end
            
            eqSpacing = 0;
        end
    end
end
    
%derivative
if which == 1
    p=input("At what point would you like to perform the derivative: ");
    
    %checking if point is in original data and spacing is even
    if ismember(p,x) && eqSpacing && find(x==p) ~= 1 && find(x==p) ~= Size
        %CDD formula
        f1 = y(find(x==p)+1);
        f2 = y(find(x==p)-1);
        der = (f1-f2)/(2*minSeg);
        
        fprintf("The derivative at this point (calculated with CDD) is %d\n", der);
    else        
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
        
        % St and Sr
        StVec=x;
        SrVec=x;

        %given formulae
        for i=1:Size
            StVec(i)=(y(i)-mean(y))^2;
            SrVec(i)=(y(i)-a0-a1*x(i)-a2*(x(i))^2-a3*(x(i))^3)^2;
        end
        
        %R^2 value
        R=(sum(StVec)-sum(SrVec))/sum(StVec);
        
        %equation and title to plot
        Y=a0+x0.*a1+(x0.^2).*a2+(x0.^3).*a3;
        Title=['Polynomial, y=' num2str(a0) '+' num2str(a1) 'x+' num2str(a2) 'x^{2}+' num2str(a3) 'x^{3}, R^{2}=' num2str(R)];
 
        figure;
        plot(x,y,'.'); %raw data
        hold on;
        plot(x0,Y); %estimated regression line
        xlabel('x');
        ylabel('y');
        %axis to fit min and max data
        axis([min(MyMatrix(:,1)) max(MyMatrix(:,1)) min(MyMatrix(:,2)) max(MyMatrix(:,2))]);
        %predetermined title
        title(Title);
        
        %CDD formula
        f1 = a0+(p+minSeg)*a1+((p+minSeg)^2)*a2+((p+minSeg)^3)*a3;
        f2 = a0+(p-minSeg)*a1+((p-minSeg)^2)*a2+((p-minSeg)^3)*a3;
        der = (f1-f2)/(2*minSeg);
        
        fprintf("The derivative at this point (calculated with CDD) is %d\n", der);
    end
else
    p1=input("Enter first integration limit: ");
    p2=input("Enter second integration limit (note that it must be larger than the first): ");
    n=input("Enter the number of segments: ");
    
    %ensuring valid bounds
    if p1 >= p2 || p1 < x(1) || p2 > x(Size)
        disp("Invalid range");
    else
        %(b-a)/n
        h = (p2-p1)/n;
        %new data set
        X = zeros(n+1);
        X(1) = p1;
        
        for i=1:n
            X(i+1) = p1+i*h;
        end
        
        subset = 1;
        
        %checking if original data points match new set
        for i=1:(n+1)
            if ismember(x,X(i)) == 0
                subset = 0;
            end
        end
        
        if subset == 1
            sum = 0;
            for i=1:(n-1)
                sum = sum+y(find(x==(p1+i*h)));
            end
            
            f1 = y(find(x==p1));
            f2 = y(find(x==p2));
            
            %trapezoid rule formula
            int = (h/2)*(f1+f2+2*sum);
            
            fprintf("The integral between these bounds (calculated with trapezoid rule) is %d\n", int);
        else
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

            % St and Sr
            StVec=x;
            SrVec=x;

            %given formulae
            for i=1:Size
                StVec(i)=(y(i)-mean(y))^2;
                SrVec(i)=(y(i)-a0-a1*x(i)-a2*(x(i))^2-a3*(x(i))^3)^2;
            end

            %R^2 value
            R=(sum(StVec)-sum(SrVec))/sum(StVec);

            %equation and title to plot
            Y=a0+x0.*a1+(x0.^2).*a2+(x0.^3).*a3;
            Title=['Polynomial, y=' num2str(a0) '+' num2str(a1) 'x+' num2str(a2) 'x^{2}+' num2str(a3) 'x^{3}, R^{2}=' num2str(R)];

            figure;
            plot(x,y,'.'); %raw data
            hold on;
            plot(x0,Y); %estimated regression line
            xlabel('x');
            ylabel('y');
            %axis to fit min and max data
            axis([min(MyMatrix(:,1)) max(MyMatrix(:,1)) min(MyMatrix(:,2)) max(MyMatrix(:,2))]);
            %predetermined title
            title(Title);
            
            sum = 0;
            for i=1:(n-1)
                sum = sum + a0+(p1+i*h)*a1+((p1+i*h)^2)*a2+((p1+i*h)^3)*a3;
            end
            
            f1 = a0+(p1)*a1+((p1)^2)*a2+((p1)^3)*a3;
            f2 = a0+(p2)*a1+((p2)^2)*a2+((p2)^3)*a3;
            
            %trapezoid rule formula
            int = (h/2)*(f1+f2+2*sum);
            
            fprintf("The integral between these bounds (calculated with trapezoid rule) is %d\n", int);
        end
    end
end