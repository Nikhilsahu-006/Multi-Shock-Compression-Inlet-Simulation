# Multi-Shock-Compression-Inlet-Simulation
Supersonic air-breathing propulsion systems — such as those in ramjets, scramjets, and high-speed aircraft — rely on shock waves to compress incoming air prior to combustion.
This project provides a numerical simulation of the complete inlet compression process, including:

one or more oblique shock forms that are reflected shocks and terminated by oblique shock when theta max limit out of domain

A terminal normal shock to transition the flow from supersonic to subsonic

A final isentropic diffuser for additional pressure recovery

Detailed tracking of static and stagnation properties through every stage

The tool is designed to be simple, modular, and sufficiently accurate for engineering analysis or academic applications.

Problem statement:

Flow is compressed through a series of oblique shocks terminated by
a normal shock (When exceeds local theta max). M1 = 3 at STP. It is required to limit M2 below
0.1 so added a diffuser at the end. Check theta1, theta2 and compute the variation of P, T, P0 along the centerline.

# Requirement

1. python 
2. numpy , math ,  scipy ,  matplotlib .


