
% box is 200nm by 100nm, 1 tick represents 1nm
BOX_WIDTH = 200;
BOX_HEIGHT = 100;
SPACE_STEP=1;

% time step are 1/1000 s, there are N_TIME_STEPS of them in the simulation
TIME_STEP = 0.1;
N_TIME_STEPS= 1000;

% rest mass in kg
ELEC_RESS_MASS = 9.10938356 * 10^(-4);
ELEC_EFF_MASS = 0.26*ELEC_RESS_MASS;

% number of particles
NUM_PARTICLES = 10000;

% temperature in kelvin
TEMPERATURE = 300;

%Boltzmann constant in J/K
kB=1.38064852*10^(-23);

%thermal velocity using the rms speed, in units of SPACE_STEP/TIME_STEP
vth = sqrt(3*kB*TEMPERATURE/ELEC_EFF_MASS) * 10^9 * TIME_STEP/ SPACE_STEP;

% generate the initial list of particles in the form of
% Px1 Py1 Vx1 Vy1
% Px2 Py2 Vx2 Vy2 and so on for a total of NUM_PARTICLES rows
particles=vec2mat(cell2mat(arrayfun(@(r)[BOX_WIDTH*rand, ...
    BOX_HEIGHT*rand, cos(2*pi*r)*vth, sin(2*pi*r)*vth],...
    rand(1,NUM_PARTICLES),'UniformOutput',false)),4);


for t_ind = 2:N_TIME_STEPS
     particles(:,:,t_ind) = vec2mat(cell2mat(arrayfun(@(row) ...
         updatePhaseSpace(particles(row,:,t_ind-1)),1:NUM_PARTICLES,...
         'UniformOutput',false)),4);
end

xVals = permute(particles(:,1,:),[1,3,2]);
yVals = permute(particles(:,2,:),[1,3,2]);

for n_ind = 1:NUM_PARTICLES
    
    breakPoints = [];
    currentXVals=xVals(n_ind,:);
    currentyVals=yVals(n_ind,:);
    
    for t_ind = 2:N_TIME_STEPS
        if abs(currentXVals(t_ind) - currentXVals(t_ind-1))> 25
            breakPoints = [breakPoints,t_ind-1];
        end
    end
    
   breakPoints= fliplr(breakPoints);
   
   for b = breakPoints
       currentXVals=[currentXVals(1:b),NaN,NaN,currentXVals(b+1:end)];
       currentyVals=[currentyVals(1:b),NaN,NaN,currentyVals(b+1:end)];
   end

    plot(currentXVals,currentyVals)
    hold on
end
hold off
    
    
    