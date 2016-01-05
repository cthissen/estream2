# estream2: Evenly spaced streamlines (2D) 

What is it?
----------------- 
estream2 calculates evenly spaced streamlines. 

The details of the method can be found in  B. Jobard and W. Lefer, Creating Evenly-Spaced Streamlines
of Arbitrary Density, Proc. Eighth Eurographics Workshop on Visualization in
Scientific Computing, pp. 45-55, 1997.

For comments, questions, or suggestions, please email cthissen@gmail.com or 
leave a comment under the issues tab at github.com/cthissen/estream2

Christopher J. Thissen, Yale University  

Why?
-----------------
Evenly spaced streamlines better illustrate the structure of the flow field. Compare even and unevely spaced streamlines:  
<img src="https://github.com/cthissen/estream2/blob/master/fig_estream2.png" alt="alt text" width="400px" height="400px">
<img src="https://github.com/cthissen/estream2/blob/master/fig_stream2.png" alt="alt text" width="400px" height="400px">


Requirements & Installation
------------------ 
No toolboxes are required to run the program, and no installation is necessary.

Usage
------------------ 
estream2 can be run with the same input arguments as stream2. For example, with estream2.m on the current path (e.g. cd to the directory containing estream2.m), run the following lines
````
load wind
sx = 80; % set initial seed
sy = 35;
stepsize = 0.05; 
maxIter = 10000; % max iterations for a single streamline
dSep = 3; % separation distance between streamlines
XY = estream2(x(:,:,5),y(:,:,5),u(:,:,5),v(:,:,5),sx,sy,[stepsize,maxIter,dSep],'plot');
streamline(XY);
````
Additional details can be found by typing help estream2 at the command prompt.


The Latest Version
------------------ 
Details of the latest version can be found on the github project page under 
  server project page under https://github.com/cthissen/Drex-MATLAB

Contributors
------------------ 
Christopher Thissen, Yale University. cthissen@gmail.com


Feedback
------------------ 
Your comments are welcome! If you find any bugs or have feature requests report them to
Christopher Thissen, christopher.thissen@yale.edu. 

Issues can also be reported online: https://github.com/cthissen/estream2/issues


License
------------------ 
See License file.
