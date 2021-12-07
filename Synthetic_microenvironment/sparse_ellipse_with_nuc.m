function [Cell,Nucleus]=sparse_ellipse_with_nuc(ra,rb,ang,x0,y0,sparsity,rat_nuc_cell,polarity)
% Ellipse adds ellipses to the current plot
%
% ELLIPSE(ra,rb,ang,x0,y0) adds an ellipse with semimajor axis of ra,
% a semimajor axis of radius rb, a semimajor axis of ang, centered at
% the point x0,y0.
%
% The length of ra, rb, and ang should be the same. 
% If ra is a vector of length L and x0,y0 scalars, L ellipses
% are added at point x0,y0.
% If ra is a scalar and x0,y0 vectors of length M, M ellipse are with the same 
% radii are added at the points x0,y0.
% If ra, x0, y0 are vectors of the same length L=M, M ellipses are added.
% If ra is a vector of length L and x0, y0 are  vectors of length
% M~=L, L*M ellipses are added, at each point x0,y0, L ellipses of radius ra.
%
% ELLIPSE(ra,rb,ang,x0,y0,C)
% adds ellipses of color C. C may be a string ('r','b',...) or the RGB value. 
% If no color is specified, it makes automatic use of the colors specified by 
% the axes ColorOrder property. For several circles C may be a vector.
%
% ELLIPSE(ra,rb,ang,x0,y0,C,Nb), Nb specifies the number of points
% used to draw the ellipse. The default value is 300. Nb may be used
% for each ellipse individually.
%
% h=ELLIPSE(...) returns the handles to the ellipses.
%
% usage exmple: the following produces a red ellipse centered at 1,1
% and tipped down at a 45 deg axis from the x axis
% ellipse(1,2,pi/4,1,1,'r')
%
% note that if ra=rb, ELLIPSE plots a circle
%
% written by D.G. Long, Brigham Young University, based on the
% CIRCLES.m original 
% written by Peter Blattner, Institute of Microtechnology, University of 
% Neuchatel, Switzerland, blattner@imt.unine.ch
% Check the number of input arguments 
if nargin<1,
  ra=[];
end;
if nargin<2,
  rb=[];
end;
if nargin<3,
  ang=[];
end;
if nargin<5,
  x0=[];
  y0=[];
end;

C=[];
Nb=[];
% set up the default values
if isempty(ra),ra=1;end;
if isempty(rb),rb=1;end;
if isempty(ang),ang=0;end;
if isempty(x0),x0=0;end;
if isempty(y0),y0=0;end;
if isempty(Nb),Nb=20;end;
if isempty(C),C=get(gca,'colororder');end;
% work on the variable sizes
% x0=ra+1;
% y0=ra+1;
ra=ra(:);
rb=rb(:);
ang=ang(:);
Nb=Nb(:);
if isstr(C),C=C(:);end;
if length(ra)~=length(rb),
  error('length(ra)~=length(rb)');
end;
if length(x0)~=length(y0),
  error('length(x0)~=length(y0)');
end;
% how many inscribed elllipses are plotted
if length(ra)~=length(x0)
  maxk=length(ra)*length(x0);
else
  maxk=length(ra);
end;
% drawing loop
for k=1:maxk
  
  if length(x0)==1
    xpos=x0;
    ypos=y0;
    radm=ra(k);
    radn=rb(k);
    if length(ang)==1
      an=ang;
    else
      an=ang(k);
    end;
  elseif length(ra)==1
    xpos=x0(k);
    ypos=y0(k);
    radm=ra;
    radn=rb;
    an=ang;
  elseif length(x0)==length(ra)
    xpos=x0(k);
    ypos=y0(k);
    radm=ra(k);
    radn=rb(k);
    an=ang(k)
  else
    rada=ra(fix((k-1)/size(x0,1))+1);
    radb=rb(fix((k-1)/size(x0,1))+1);
    an=ang(fix((k-1)/size(x0,1))+1);
    xpos=x0(rem(k-1,size(x0,1))+1);
    ypos=y0(rem(k-1,size(y0,1))+1);
  end;
  co=cos(an);
  si=sin(an);
  the=linspace(0,2*pi,Nb(rem(k-1,size(Nb,1))+1,:)+1);

% Create cell 
x=radm*cos(the)*co-si*radn*sin(the)+xpos;
y=radm*cos(the)*si+co*radn*sin(the)+ypos;

% Create nucleus
x_nuc=radm*rat_nuc_cell*cos(the)*co-si*radn*rat_nuc_cell*sin(the)+xpos;
y_nuc=radm*rat_nuc_cell*cos(the)*si+co*rat_nuc_cell*radn*sin(the)+ypos;
x_nuc = x_nuc + rand(1,Nb+1).*x0.*0.05;
y_nuc = y_nuc + rand(1,Nb+1).*y0.*0.05;


% Add sparseNoise 
x = x + full(sprandn(1,Nb+1,sparsity/2)).*x0./7 + rand(1,Nb+1).*x0./10;
y = y + full(sprandn(1,Nb+1,sparsity/2)).*y0./7 + rand(1,Nb+1).*y0./10;

% Interpolation
mult=10;
x = [x,x,x]; y = [y,y,y];
x = interpn(1:(Nb+1)*3,x,1:1/mult:(Nb+1)*3,'cubic')';
y = interpn(1:(Nb+1)*3,y,1:1/mult:(Nb+1)*3,'cubic')';
x = x((Nb+1)*mult:2*(Nb+1)*mult); y = y((Nb+1)*mult:2*(Nb+1)*mult); 

% Interpolation
mult=10;
x_nuc = [x_nuc,x_nuc,x_nuc]; y_nuc = [y_nuc,y_nuc,y_nuc];
x_nuc = interpn(1:(Nb+1)*3,x_nuc,1:1/mult:(Nb+1)*3,'cubic')';
y_nuc = interpn(1:(Nb+1)*3,y_nuc,1:1/mult:(Nb+1)*3,'cubic')';
x_nuc = x_nuc((Nb+1)*mult:2*(Nb+1)*mult); y_nuc = y_nuc((Nb+1)*mult:2*(Nb+1)*mult); 

% Apply polarity to the nucleus
x_pol_mov = co*radm*polarity/4 + si*radn*polarity/4;
y_pol_mov = si*radm*polarity/4 + co*radn*polarity/4;
r = rand;
if r<=0.5
    x_nuc = x_nuc+x_pol_mov; y_nuc = y_nuc+y_pol_mov;
else
    x_nuc = x_nuc-x_pol_mov; y_nuc = y_nuc-y_pol_mov;
end

% Check values
min_x = min(x);min_y = min(y);
x = x-min_x+1; y = y-min_y+1;
x(isnan(x)) = 1; x(isinf(x)) = 1;
y(isnan(y)) = 1; y(isinf(y)) = 1;
x_nuc = x_nuc-min_x+1; y_nuc = y_nuc-min_y+1;
x_nuc(isnan(x_nuc)) = 1; x_nuc(isinf(x_nuc)) = 1;
y_nuc(isnan(y_nuc)) = 1; y_nuc(isinf(y_nuc)) = 1;

% Create cell
Cell = poly2mask(x,y,ceil(max(y)/2)*2,ceil(max(x)/2)*2);
Nucleus = poly2mask(x_nuc,y_nuc,ceil(max(y)/2)*2,ceil(max(x)/2)*2);

if all(Cell==0)
    Cell(1,1) = 1;
end
if all(Nucleus==0)
    Nucleus(1,1) = 1;
end

Cell(Nucleus>0) = 1;
% figure; imshow(label2rgb(Cell+Nucleus))

end
