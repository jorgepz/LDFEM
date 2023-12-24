//ms = 0.05 ;
ms = 1.0 ;

Point(1) = {0.0,0.0,0.0 , ms }; //
Point(2) = {1.0,0.0,0.0 , ms }; //
Point(3) = {1.0,1.0,0.0 , ms }; //
Point(4) = {0.0,1.0,0.0 , ms }; //
Point(5) = {0.0,0.0,1.0 , ms }; //
Point(6) = {1.0,0.0,1.0 , ms }; //
Point(7) = {1.0,1.0,1.0 , ms }; //
Point(8) = {0.0,1.0,1.0 , ms }; //

Line(1)  = {1,2} ;
Line(2)  = {2,3} ;
Line(3)  = {3,4} ;
Line(4)  = {4,1} ;
Line(5)  = {5,6} ;
Line(6)  = {6,7} ;
Line(7)  = {7,8} ;
Line(8)  = {8,5} ;
Line(9)   = {1,5} ;
Line(10)  = {2,6} ;
Line(11)  = {3,7} ;
Line(12)  = {4,8} ;

// z=0
Line Loop(1) = {1,2,3,4} ;
// z=1
Line Loop(2) = {-5,-8,-7,-6} ;
// y=0
Line Loop(3) = {-1,9,5,-10} ;
// x=1
Line Loop(4) = {-2,10,6,-11} ;
// y=1
Line Loop(5) = {-3,11,7,-12} ;
// x=0
Line Loop(6) = {-4,12,8,-9} ;

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};

Physical Point(1) = {1,2,3,4,5,6,7,8} ;

Physical Surface(1) =  {2,5} ;
Physical Surface(2) =  {6} ;
Physical Surface(3) =  {3} ;
Physical Surface(4) =  {1} ;
Physical Surface(5) =  {4} ;

Physical Volume(1)   = {1} ;

//Field[50] = MathEval; //Generate Field
//Field[50].F = "1";
//Background Field = 50;
