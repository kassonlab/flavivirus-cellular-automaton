function [CumX, CumY] = Calculate_CDF(WaitTimeList)
%Changes WaitTimeList to proportion hemifused over time and makes cum
%dist(Not normalized)

    WaitTimeList = sort(WaitTimeList);
    
    Y = 1:length(WaitTimeList);
    New_Pt = 0;
    
    for i = 1:length(Y)
       if i == length(Y)
            New_Pt = New_Pt +1;
           CumX(New_Pt) = WaitTimeList(i);
           CumY(New_Pt) = Y(i);
           continue
       end
       if WaitTimeList(i) == WaitTimeList(i+1)
           continue
       else
           New_Pt = New_Pt +1;
           CumX(New_Pt) = WaitTimeList(i);
           CumY(New_Pt) = Y(i);
       end
    end
 
end
