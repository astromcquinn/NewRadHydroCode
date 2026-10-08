#!/usr/bin/env python3
"""Photon deficit along the grid axes from a single source (quirk B16,
docs/port_notes.md). Usage: axis_deficit.py RTDRIVER_OUTDIR... (runs of
tests/problems/uniform_1src, source cell 16 = index 15, N = 32). Prints, per
step, the most negative x_HII on the +-x, +-y, +-z axis cells minus the mean of
the shell at the same radius."""
import numpy as np, sys, glob, os
N=32; s0=15  # source cell (0-based), centre 15.5
for run in sys.argv[1:]:
    for step in (5,10,20,30):
        p=os.path.join(run,"step_%04d"%step,"fHI.bin")
        if not os.path.exists(p): continue
        x=1-np.fromfile(p,'<f4').reshape(N,N,N,order='F')
        i=np.arange(N)
        X,Y,Z=np.meshgrid(i,i,i,indexing='ij')
        r=np.sqrt((X-s0)**2+(Y-s0)**2+(Z-s0)**2)
        out=[]
        for name,mask in (("x",(Y==s0)&(Z==s0)),("y",(X==s0)&(Z==s0)),("z",(X==s0)&(Y==s0))):
            d=[]
            for rr in range(2,14):
                sh=(np.abs(r-rr)<0.5)
                on=mask&(np.round(r)==rr)
                if on.any(): d.append(x[on].mean()-x[sh&~mask].mean())
            out.append("%s:%+.3f"%(name,min(d)))
        # front radius
        R=((x>0.5).sum()*3/(4*np.pi))**(1/3)
        print("%-8s step %2d R=%.1f  min(axis - shell mean) x_HII  %s"%(os.path.basename(run),step,R,"  ".join(out)))
