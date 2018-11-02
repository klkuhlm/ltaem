Intersection of conics matlab functions
For more info: www.pigei.com or pierluigi.taddei@polimi.it

%%example usage
%%the homogeneous representation of a conic is a matrix
%% m = [A C D; C B E; D E F] that represents the equation
%% A x^2 + B y^2 + 2C xy + 2D x + 2Ey + F = 0


%a circle centered in the origin
E1 = [1 0 0; 0 1 0; 0 0 -3]

%an ellipse centered in the origin
E2 = [1 0 0; 0 3 0; 0 0 -6]

%get the four homogeneous intersections
P = intersectConics(E1, E2)

%plot the normalized points
plot(P(1,:) ./ P(3,:) , P(2,:) ./ P(3,:), 'ro');







%%%%%%%%%changes

v.1.0.2: bug fixes (characteristic polynom)
v.1.0.4: bug fixes (special case of linear equations)