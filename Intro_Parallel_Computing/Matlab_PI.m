clear all;clc
TOTROUNDS=500000000;
%%% SERIAL %%%
tic
score=0;
for dart=1:TOTROUNDS
    X2=((2*rand)-1)^2;
    Y2=((2*rand)-1)^2;
    if (X2+Y2 < 1)
        score=score+1;
    end
end
pi=4*score/TOTROUNDS
toc
%%% PARFOR %%%
tic 
score=0;
parfor dart=1:TOTROUNDS
    X2=((2*rand)-1)^2;
    Y2=((2*rand)-1)^2;
    if (X2+Y2 < 1)
        score=score+1;
    end
end
pi=4*score/TOTROUNDS
toc
%%% SPMD  %%%
tic
score=0;
spmd
    ROUNDS=int64(TOTROUNDS/numlabs);
    for dart=1:ROUNDS
     X2=((2*rand)-1)^2;
     Y2=((2*rand)-1)^2;
     if (X2+Y2 < 1)
        score=score+1;
     end
    end
    if (labindex ~= 1 )
      labSend(score,1) 
    end
 
   if (labindex == 1)
      for ww=2:numlabs
          
         score=score+labReceive(ww);
      end
   end
   labBarrier
   if (labindex == 1)
      pi=4*score/(ROUNDS*numlabs) 
   end
end
toc
