clf;
global n; 
n=75; 
epsilon = 0.001;
global fuzz;
fuzz = 0.5;
global numberOfClisters;
numberOfClisters = 3;
global belongBorder
belongBorder = 0.4;

%some ][ynia

xyu = zeros(75,2);
xyu1 = zeros(75,2);
xyu2 = zeros(75,2);

%crete arrays of zeroes
Sorted = rand(3,n,2);

%create clasters elements

clasterElementA = [0.9,0.5];
clasterElementB = [0.1,0.1];
clasterElementC = [0.1,0.9];
%intialize data

plot(Sorted(1,:,1),Sorted(1,:,2),'g*');hold on
plot(Sorted(2,:,1),Sorted(2,:,2),'r*');hold on
plot(Sorted(3,:,1),Sorted(3,:,2),'b*');hold on

plot(clasterElementA(1),clasterElementA(2),'X');hold on
plot(clasterElementB(1),clasterElementB(2),'X');hold on
plot(clasterElementC(1),clasterElementC(2),'X');



Stop = false;
numberOFtrying = 0;

fuzzyCoeficientMatrix = zeros(numberOfClisters,n,3);
distanceToClaster = zeros(numberOfClisters,n,3);

pause(5);
while (~Stop)&& numberOFtrying < 10000
   numberOFtrying = numberOFtrying+1;
   for j =1:numberOfClisters
    for i = 1:n
      [fuzzyCoeficientMatrix(j,i,1),fuzzyCoeficientMatrix(j,i,2),fuzzyCoeficientMatrix(j,i,3),distanceToClaster(j,i,1),distanceToClaster(j,i,2),distanceToClaster(j,i,3)] = calculateFuzzyCoeficent(Sorted(j,i,:),clasterElementA,clasterElementB,clasterElementC);
    end
   end    
  
   summurize = 0;
   for j = 1:numberOfClisters
     for i = 1:n
        summurize = summurize + distanceToClaster(j,i,1)*fuzzyCoeficientMatrix(j,i,1);
        summurize = summurize + distanceToClaster(j,i,2)*fuzzyCoeficientMatrix(j,i,2);
        summurize = summurize + distanceToClaster(j,i,3)*fuzzyCoeficientMatrix(j,i,3);
      end
   end
    xyu(numberOFtrying,:) = clasterElementA;
 xyu1(numberOFtrying,:) = clasterElementB;
 xyu2(numberOFtrying,:) = clasterElementC;
   
   if(numberOFtrying > 1)       
      Stop =  (abs(oldDistance - summurize)<epsilon);
   end
   
   oldDistance = summurize;


   %MISTAKE IN N
 [clasterElementA(1),clasterElementA(2)] =  calculateClasters(fuzzyCoeficientMatrix,Sorted,1);
 [clasterElementB(1),clasterElementB(2)] =  calculateClasters(fuzzyCoeficientMatrix,Sorted,2);
 [clasterElementC(1),clasterElementC(2)] =  calculateClasters(fuzzyCoeficientMatrix,Sorted,3);


  
end

clf;

A = zeros(2*n,2);
B = zeros(2*n,2);
C = zeros(2*n,2);
countA=0;
countB=0;
countC=0;


for j=1:numberOfClisters
    for i=1:n
       ClasterAKoef = fuzzyCoeficientMatrix(j,i,1);
       ClasterBKoef = fuzzyCoeficientMatrix(j,i,2);
       ClasterCKoef = fuzzyCoeficientMatrix(j,i,3);
       if(ClasterAKoef>ClasterBKoef)&&(ClasterAKoef>ClasterCKoef)
        if ClasterAKoef > belongBorder
          countA = countA +1;
          A(countA,:) = Sorted(j,i,:);
        end
       end
       if(ClasterBKoef>ClasterAKoef)&&(ClasterBKoef>ClasterCKoef)
         if ClasterBKoef >belongBorder
          countB = countB +1;
          B(countB,:) = Sorted(j,i,:);
         end
       end      
       if(ClasterCKoef>ClasterAKoef)&&(ClasterCKoef>ClasterBKoef)
         if ClasterCKoef > belongBorder
          countC = countC +1;
          C(countC,:) = Sorted(j,i,:);
         end
       end
    end
end

A=A(any(A,2),:);
B=B(any(B,2),:);
C=C(any(C,2),:);

plot(A(:,1),A(:,2),'g*');hold on
plot(B(:,1),B(:,2),'r*');hold on
plot(C(:,1),C(:,2),'b*');hold on

plot(clasterElementA(1),clasterElementA(2),'X');hold on
plot(clasterElementB(1),clasterElementB(2),'bX');
plot(clasterElementC(1),clasterElementC(2),'bX');


function [fuzzyKoeficientA,fuzzyKoeficientB,fuzzyKoeficientC,distanceA,distanceB,distanceC]  = calculateFuzzyCoeficent(Array,clasterElementA,clasterElementB,clasterElementC)
     global fuzz;
     
     distanceA = myDistance(Array(1,1,1),Array(1,1,2),clasterElementA(1),clasterElementA(2));
     distanceB = myDistance(Array(1,1,1),Array(1,1,2),clasterElementB(1),clasterElementB(2));
     distanceC = myDistance(Array(1,1,1),Array(1,1,2),clasterElementC(1),clasterElementC(2));
     
     fuzzyKoeficientA = (1/distanceA)^(2/fuzz-1);
     fuzzyKoeficientB = (1/distanceB)^(2/fuzz-1);
     fuzzyKoeficientC = (1/distanceC)^(2/fuzz-1);
     
     sumForZnamenatel = fuzzyKoeficientA+fuzzyKoeficientB+fuzzyKoeficientC;
     
     fuzzyKoeficientA = fuzzyKoeficientA/(sumForZnamenatel);
     fuzzyKoeficientB = fuzzyKoeficientB/(sumForZnamenatel);
     fuzzyKoeficientC = fuzzyKoeficientC/(sumForZnamenatel);
     
     
end

function [X,Y] = calculateClasters(fuzzyCoeficientMatrix,array,numberOfCluxter)
  global fuzz;
  global numberOfClisters;
  global n;
  global belongBorder;
  
  clasterX = zeros(2);
  clasterY = zeros(2);
    for j = 1:numberOfClisters
      for i = 1:n
       if  fuzzyCoeficientMatrix(j,i,numberOfCluxter) > belongBorder
            clasterX(1) = clasterX(1) + ((fuzzyCoeficientMatrix(j,i,numberOfCluxter)^fuzz)*array(j,i,1));
            clasterX(2) = clasterX(2) + (fuzzyCoeficientMatrix(j,i,numberOfCluxter)^fuzz);
            clasterY(1) = clasterY(1) + ((fuzzyCoeficientMatrix(j,i,numberOfCluxter)^fuzz)*array(j,i,2));
            clasterY(2) = clasterY(2) + (fuzzyCoeficientMatrix(j,i,numberOfCluxter)^fuzz);
       end
      end
    end
    X = clasterX(1)/clasterX(2);
    Y = clasterY(1)/clasterY(2);
    
end

function y=myDistance(x,y,x1,y1)
    y = sqrt((x-x1)^2+(y-y1)^2);
end



