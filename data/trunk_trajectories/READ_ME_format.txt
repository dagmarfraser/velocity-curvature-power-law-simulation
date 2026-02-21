*** Notes on the trajectories files ***

*** 2020_03_13_Chishuru ***

10 rows of markers (1 = head base, 10 = tip)
number of markers per row  = (row:markers) = (1:6, 2:4, 3:6, 4:6, 5:6, 6:5, 7:5, 8:4, 9:3, 10:3)

*** File format for the markers positions ***

lines = time frame (frame rate = 100 Hz)
columns = x,y,z coordinates of marker i from row j 
	(row 1, marker 1, x) (row 1, marker 1, y) … (row 1, marker 2, x) … (last row, last marker, z)
marker 1 is on the right side of the animal

*** File format for the ellipses parameters ***

lines = time frame (frame rate = 100 Hz)
columns = 8 parameters of each ellipse along the trunk, from 1 (head base) to 9 or 10 (tip)
	(param1_1) … (param8_1) … (param1_10) … (param8_10)
parameters given in the following order: 
	(1) large_axis (mm)
	(2) ratio of large over small axis
	(3) alpha (radians)
	(4) beta (radians)
	(5) gamma (radians)
	(6), (7), (8) x,y,z coordinates (mm) of ellipse center
the angles alpha, beta and gamma give the orientation of the ellipse with the following rotation matrix:
	       Rz = [cos(gamma), -sin(gamma), 0; ...
            	   sin(gamma), cos(gamma), 0; ... 
          		    0, 0, 1];
      
       	       Ry = [cos(beta), 0, sin(beta);  ...
          		    0, 1, 0;  ...
          		    -sin(beta), 0, cos(beta)];

    	       Rx = [1, 0, 0;  ...
           		    0, cos(alpha), -sin(alpha);  ...
            	    0, sin(alpha), cos(alpha)];

     	       R = Rz*Ry*Rx;

	       So that the ellipse principal axes are R*[a;0;0] and R*[0;b;0] 