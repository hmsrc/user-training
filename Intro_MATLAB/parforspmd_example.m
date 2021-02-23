clear all;clc
parpool('local',4)

ARRAY=rand(100,1);
%%
tic
spmd
    
   %disp('I am worker')
   %labindex
   LOW=(labindex-1)*length(ARRAY)/numlabs+1;
   HIGH=length(ARRAY)/numlabs*labindex;
   tot=sum(ARRAY((LOW:HIGH)));    
   if (labindex ~= 1 )
      labSend(tot,1) 
   end
 
   if (labindex == 1)
      for ww=2:numlabs
          
         tot=tot+labReceive(ww);
      end
   end
   labBarrier
   if (labindex == 1)
    %  disp('big total') 
      tot 
   end
end
toc
%%
tic
tot2=0; 
 parfor i=1:length(ARRAY)
     
    tot2=tot2+ARRAY(i); 
   
  
 end
tot2
toc
%%
tic
tot3=0;
for  i=1:length(ARRAY)
     
    tot3=tot3+ARRAY(i); 
     
  
 end
tot3
toc
%%
tic
sum(ARRAY(:))
toc

%%

tic;parfor i=1:10;pause(1);end;toc
%%
tic;parfor i=1:12;pause(1);end;toc

