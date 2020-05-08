/* particles.c Simulation of an exploding bag of bouncy balls / particles

Particles explode outward, fall under the weight of gravity, and bounce off the floor.

Uses simple kinematic definitions from introductory physics:
    velocity = dx / dt
    acceleration = dv / dt

Uses the known acceleration from gravity:
    g = -9.8 m/s^2 (in the y direction)
    
And takes advantage of the physics of an "inelastic collision" between a bouncy ball and the floor.
    x-velocity is unchanged (velx = velx)
    y-velocity reverses direction but keeps magnitude (vely = -vely)
    
Basically just repeats these simple equations many, many times.

Compile with:
    gcc -Wall -fopenmp particles.c -o particles.o -lm

References:
    * C Arrays tutorial: https://beginnersbook.com/2014/01/c-arrays-example/
    * How to do multiple outputs in C: https://www.geeksforgeeks.org/how-to-return-multiple-values-from-a-function-in-c-or-cpp/
    * Random numbers in C https://c-for-dummies.com/blog/?p=1458

Created by Scott Feister on April 22, 2020
*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>

#define NPARTS 1000 // Number of particles
#define NSTEPS 600 // Number of steps to record (60 fps for 10 seconds playback = 600 steps)
#define NSUBSTEPS 5000 // Number of substeps per step (increases physical accuracy; also increases runtime)

void new_particle(double* x, double* y, double* velx, double* vely);

void new_particle(double* x, double* y, double* velx, double* vely) {
    // Randomly generates a new particle; that is, assigns random values for x position and velocity
    // Note that the expression '(float)rand()/RAND_MAX' gives a random number between 0 and 1.
    *x = 5.0 * (float)rand()/RAND_MAX; // Initial x value for particles will be between x = 0 meter and x = +5.0 meters
    *y = 10.0 + 5.0 * (float)rand()/RAND_MAX; // Initial y value for particles will be between y = +10.0 meters and y = +15.0 meters
    *velx = -10.0 + 20.0 * (float)rand()/RAND_MAX; // initial x-velocity will be between -10.0 and +10.0 meters/second
    *vely = 20.0 * (float)rand()/RAND_MAX; // initial y-velocity will be between 0 and +20.0 meters/second
}

int main() {
    printf("Assigning / allocating variables and arrays.\n");

    // Set final time
    const double tfinal = 10.0; // We will simulate 10 seconds of particle time
    
    // Initialize variables relating to each individual particle
    double x; // x-value of particle
    double y; // y-value of particle
    double velx; // x-component particle velocity
    double vely; // y-component particle velocity
    const double accx = 0.0; // x-component particle acceleration: none, since gravity is in the y direction.
    const double accy = -9.8; // y-component particle acceleration: for gravity on Earth, always -9.8 meters/second^2

    // Allocate particle arrays to store results
    float xarr[NPARTS*NSTEPS]; // array of x-values of particles, at each recorded step
    float yarr[NPARTS*NSTEPS]; // array of y-values of particles, at each recorded step
    
    // Compute timestep
    const double dt = tfinal / NSTEPS / NSUBSTEPS; // time step of simulation

    // Begin simulation
    printf("Running particle simulation...\n");
    double t0 = omp_get_wtime(); // Start timer
    int i, j, step;
    
    // Iterate over particles
    for (i=0; i<NPARTS;i++)
    {
        // Make particle
        new_particle(&x, &y, &velx, &vely);
       
        // Iterate over timesteps, for this particle
        for (step=0; step<NSTEPS;step++)
        {                                
            // Save position values into array
            xarr[i*NSTEPS + step] = x;
            yarr[i*NSTEPS + step] = y;

            for (j=0; j<NSUBSTEPS;j++)
            {
                // If particle hits the ground (y=0), bounce back up (reverse y-velocity)
                if (y < 0 && vely < 0)
                {
                    vely = -vely;
                }
                
                // Advance timestep according to simple Newtonian physics equations: v=dx/dt and a=dv/dt.
                x += dt * velx; // Compute particle's x-position for next step
                y += dt * vely; // Compute particle's y-position for next step
                velx += dt * accx; // Compute particle's x-velocity for next step
                vely += dt * accy; // Compute particle's y-velocity for next step
            }
        }
    }
    
    double t1 = omp_get_wtime(); // Stop timer
    printf("Particle simulation complete!\n");
    printf("Simulation execution time: %f secs\n", t1 - t0);
    
    // Store x and y position arrays to CSV files 'xarr.txt' and 'yarr.txt'
    //   Each row will be a particle
    //   Each column will be a step
    FILE *fpx, *fpy;
    fpx = fopen("xarr.txt","w");
    fpy = fopen("yarr.txt","w");
    for (i=0; i<NPARTS;i++)
    {
        for (step=0; step<NSTEPS;step++)
        {
            fprintf(fpx, "%f,", xarr[i*NSTEPS + step]);
            fprintf(fpy, "%f,", yarr[i*NSTEPS + step]);
        }
        fprintf(fpx,"\n");
        fprintf(fpy,"\n");
    }
    fclose(fpx);
    fclose(fpy);
    
    printf("Outputs saved to 'xarr.txt' and 'yarr.txt'.\n");
    printf("Visualize results by calling 'python particles.py'.\n");
    return 0;
}
