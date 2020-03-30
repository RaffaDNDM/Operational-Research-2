\ENCODING=ISO-8859-1
\Problem name: TSP

Minimize
 obj1: 4727 x1,2 + 1204 x1,3 + 6363 x1,4 + 3657 x1,5 + 4725 x2,1 + 3586 x2,3
       + 2011 x2,4 + 1840 x2,5 + 1203 x3,1 + 3587 x3,2 + 5162 x3,4 + 2458 x3,5
       + 6361 x4,1 + 2011 x4,2 + 5161 x4,3 + 2797 x4,5 + 3656 x5,1 + 1841 x5,2
       + 2457 x5,3 + 2799 x5,4 + 0 u1
Subject To
 Degreemin_constraint(1)#0:  x1,1 + x2,1 + x3,1 + x4,1 + x5,1  = 1
 Degreemin_constraint(2)#1:  x1,2 + x2,2 + x3,2 + x4,2 + x5,2  = 1
 Degreemin_constraint(3)#2:  x1,3 + x2,3 + x3,3 + x4,3 + x5,3  = 1
 Degreemin_constraint(4)#3:  x1,4 + x2,4 + x3,4 + x4,4 + x5,4  = 1
 Degreemin_constraint(5)#4:  x1,5 + x2,5 + x3,5 + x4,5 + x5,5  = 1
 Degreemout_constraint(1)#5: x1,1 + x1,2 + x1,3 + x1,4 + x1,5  = 1
 Degreemout_constraint(2)#6: x2,1 + x2,2 + x2,3 + x2,4 + x2,5  = 1
 Degreemout_constraint(3)#7: x3,1 + x3,2 + x3,3 + x3,4 + x3,5  = 1
 Degreemout_constraint(4)#8: x4,1 + x4,2 + x4,3 + x4,4 + x4,5  = 1
 Degreemout_constraint(5)#9: x5,1 + x5,2 + x5,3 + x5,4 + x5,5  = 1
 BigmM_constraint(1)#10:     5 x2,2 <= 4
 BigmM_constraint(2)#11:     5 x2,3 + u2 - u3 <= 4
 BigmM_constraint(3)#12:     5 x2,4 + u2 - u4 <= 4
 BigmM_constraint(4)#13:     5 x2,5 + u2 - u5 <= 4
 BigmM_constraint(5)#14:     5 x3,2 - u2 + u3 <= 4
 BigmM_constraint(6)#15:     5 x3,3 <= 4
 BigmM_constraint(7)#16:     5 x3,4 + u3 - u4 <= 4
 BigmM_constraint(8)#17:     5 x3,5 + u3 - u5 <= 4
 BigmM_constraint(9)#18:     5 x4,2 - u2 + u4 <= 4
 BigmM_constraint(10)#19:    5 x4,3 - u3 + u4 <= 4
 BigmM_constraint(11)#20:    5 x4,4 <= 4
 BigmM_constraint(12)#21:    5 x4,5 + u4 - u5 <= 4
 BigmM_constraint(13)#22:    5 x5,2 - u2 + u5 <= 4
 BigmM_constraint(14)#23:    5 x5,3 - u3 + u5 <= 4
 BigmM_constraint(15)#24:    5 x5,4 - u4 + u5 <= 4
 BigmM_constraint(16)#25:    5 x5,5 <= 4
Bounds
      x1,1 = 0
 0 <= x1,2 <= 1
 0 <= x1,3 <= 1
 0 <= x1,4 <= 1
 0 <= x1,5 <= 1
 0 <= x2,1 <= 1
      x2,2 = 0
 0 <= x2,3 <= 1
 0 <= x2,4 <= 1
 0 <= x2,5 <= 1
 0 <= x3,1 <= 1
 0 <= x3,2 <= 1
      x3,3 = 0
 0 <= x3,4 <= 1
 0 <= x3,5 <= 1
 0 <= x4,1 <= 1
 0 <= x4,2 <= 1
 0 <= x4,3 <= 1
      x4,4 = 0
 0 <= x4,5 <= 1
 0 <= x5,1 <= 1
 0 <= x5,2 <= 1
 0 <= x5,3 <= 1
 0 <= x5,4 <= 1
      x5,5 = 0
      u1 = 0
 0 <= u2 <= 3
 0 <= u3 <= 3
 0 <= u4 <= 3
 0 <= u5 <= 3
Binaries
 x1,1  x1,2  x1,3  x1,4  x1,5  x2,1  x2,2  x2,3  x2,4  x2,5  x3,1  x3,2  x3,3 
 x3,4  x3,5  x4,1  x4,2  x4,3  x4,4  x4,5  x5,1  x5,2  x5,3  x5,4  x5,5 
Generals
 u1  u2  u3  u4  u5 
End
