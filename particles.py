#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
particles.py: Make a movie visualization of the bouncy particles simulation

Generates a .mp4 movie file from your simulation raw outputs 'xarr.txt' and 'yarr.txt'.

Before running this script, first compile and run 'particles.c'

Created by Scott Feister on Mon Apr 13 16:01:15 2020
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def update_scat(step, data, scat):
    scat.set_offsets(data[:,:,step].T)
    return scat,

if __name__ == "__main__":
    xarr = np.genfromtxt("xarr.txt", delimiter=",")
    yarr = np.genfromtxt("yarr.txt", delimiter=",")  
    data = np.stack([xarr, yarr])
    nparts = xarr.shape[0]
    nsteps = xarr.shape[1]
    rgb = np.random.rand(nparts,3) # Randomize particle color
    
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=60, metadata=dict(artist='Me'), bitrate=8400)

    fig, ax = plt.subplots()
    
    scat = ax.scatter(xarr[:,0], yarr[:,0], s=1, c=rgb)
    ax.axhline(y = 0) # ground indicator
    ax.set_xlim(-30, 50)
    ax.set_ylim(-1, 40)
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_title('Bouncy Particles Movie')
    ax.set_aspect('equal')
    scat_ani = animation.FuncAnimation(fig, update_scat, nsteps, fargs=(data, scat), blit=True)
    
    print("Animating movie (can take a minute or two)...")
    scat_ani.save('particles.mp4', writer=writer)
    print("Simulation movie saved to 'particles.mp4'")
    print("Copy the movie file to your laptop to watch it!")
