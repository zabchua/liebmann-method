// Laplace Equation using Liebmann Over Relaxation Method
// Chua, Mary Elizabeth E. Chua 

// INSERT DIMENSIONS HERE
rows = 70
columns = 70

u = zeros(columns+2,rows+2);
uold = zeros(columns+2,rows+2);

ea = zeros(columns+2,rows+2);
es = 1;

// INSERT THE TEMPERATURE OF THE SIDES
u(columns+2,:) = 50;
u(:,rows+2) = 100;
u(1,:) = 75;
u(:,1) = 0;

// LAMBDA
lambda = 1.5;

iter = 1;
check = 0;

function B = rot90(A,k)
    if k == 1
        A = A.';
        B = A(rows+2:-1:1,:);
    elseif k == 2
        B = A(columns+2:-1:1,rows+2:-1:1);
    elseif k == 3
        B = A(columns+2:-1:1,:);
        B = B.';
    else
        B = A;
    end
endfunction

while (iter == 1 || check == 0)
    check =1;
 
    for (j=2:rows+1)
        for (i=2:columns+1)
            
            temp = (u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1))/4;
//            mprintf("at $T_{%d%d}$,\n $$T_{%d%d} = \\dfrac{%f + %f + % f + %f}{4} = %f $$\n", i-1, j-1,i-1, j-1, u(i+1, j), u(i-1, j), u(i, j+1), u(i, j-1), temp)
            u(i,j) = lambda*temp + (1-lambda)*uold(i,j);
//            mprintf("applying the over-relaxation,\n")
//            mprintf("$$T_{%d%d} = %0.1f(%f) + (1-%0.1f)%f = \\fbox{T_{%d%d}= %f} $$\n", i-1, j-1, lambda, temp, lambda, uold(i,j), i-1, j-1, u(i,j))
            
            ea(i,j) = abs((u(i,j)-uold(i,j))/u(i,j))*100;
//            mprintf("computing for the $\\epsilon_a$,\n")
//            mprintf("$$\\epsilon_{a(%d,%d)}=\\left| \\dfrac{%f - %f}{%f} \\right| 100 = %f$$\n", i-1, j-1, u(i,j), uold(i,j), u(i,j), ea(i,j))
            if(ea(i,j)>es)then
                check=0;
            end
            
            uold(i,j) = u(i,j);
//            mprintf("\n")
        end
//        mprintf("\n")
    end
    
//    for (i=2:columns+1)
//        for (j=rows+1:-1:2)
//            mprintf("$$T_{%d%d}= %f$$\n", i-1, j-1, u(i,j))
//        end
//    end
//    
//    for (i=2:columns+1)
//        for (j=rows+1:-1:2)
//            mprintf("$$\\epsilon_{a%d%d}= %f$$\n", i-1, j-1, ea(i,j))
//        end
//    end    
    if(iter<=30) then
        mprintf("AT ITERATION %d\n", iter)
        format(8)
        disp(rot90(u, 1))
//      disp(ea)
        mprintf("\n\n\n")
    end
    iter = iter + 1;
end
disp(iter)










