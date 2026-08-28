#!/usr/bin/env python3
"""Regenerate the PNS documentation figures from the installed scoring code.

From the source root:
    PYTHONPATH=src python examples/plot_pns_kernels.py
"""

from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

from nucleosuite.scoring.pns import precompute_distributions


def save(fig, directory, stem):
    for ext in ("png", "svg"):
        fig.savefig(directory / f"{stem}.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def geometry_figure(directory, mode=167):
    lengths = [120, 167, 180]
    colors = ["#2b9b74", "#cf3e4e", "#377bb5"]
    signed, positive = precompute_distributions(lengths, mode)
    fig, axes = plt.subplots(3, 3, figsize=(11.5, 7.2), sharex=True,
                             gridspec_kw={"height_ratios": [0.75, 2, 2]}, layout="constrained")
    for col, (length, color) in enumerate(zip(lengths, colors)):
        start = -max(0, mode-length)
        x = start + np.arange(len(signed[length]))
        axes[0, col].plot([start, x[-1]], [0.28, 0.28], color="0.65", lw=2)
        axes[0, col].plot([0, length-1], [0.65, 0.65], color=color, lw=6, solid_capstyle="butt")
        axes[0, col].text((length-1)/2, 0.89, f"{length} bp fragment", ha="center", color=color, weight="bold")
        axes[0, col].text((length-1)/2, 0.02, f"{len(x)} bp scoring support", ha="center", color="0.35", fontsize=10)
        axes[0, col].set_ylim(-0.1, 1.15)
        axes[0, col].axis("off")
        for row, values in [(1, signed[length]), (2, positive[length])]:
            ax=axes[row,col]
            ax.axvspan(0, length-1, color="0.94", zorder=0)
            ax.axvline(0, color="0.5", lw=0.8, ls="--")
            ax.axvline(length-1, color="0.5", lw=0.8, ls="--")
            ax.axhline(0, color="0.35", lw=0.8)
            ax.plot(x,values,color=color,lw=2.0)
            ax.set_xlim(-55,235)
            ax.xaxis.set_major_locator(MultipleLocator(50))
            ax.xaxis.set_minor_locator(MultipleLocator(10))
            ax.spines[['top','right']].set_visible(False)
            ax.grid(axis='y',alpha=0.2)
        axes[1,col].set_ylim(-2.15,2.15)
        axes[2,col].set_ylim(-0.15,4.15)
        axes[1,col].text(0.97,0.96,"+100 / −100 mass",ha="right",va="top",fontsize=9,transform=axes[1,col].transAxes)
        axes[2,col].set_xlabel("Position from fragment start (bp)")
    axes[1,0].set_ylabel("PNS contribution per base")
    axes[2,0].set_ylabel("posPNS contribution per base")
    fig.suptitle("PNS fragment geometry and native kernels · protected-DNA mode 167 bp",fontsize=14,weight="bold")
    save(fig,directory,"pns_kernels_120_167_180_mode167")


def adaptation_figure(directory, mode=167):
    lengths=[137,152,167,182,197]
    colors=["#377bb5","#9467bd","#cf3e4e","#e29331","#2b9b74"]
    signed,_=precompute_distributions(lengths,mode)
    fig=plt.figure(figsize=(11.5,7.0),layout="constrained")
    gs=fig.add_gridspec(3,2,width_ratios=[1.45,1],height_ratios=[0.85,1,1])
    spans=fig.add_subplot(gs[0,0]); wave=fig.add_subplot(gs[1:,0]); amplitude=fig.add_subplot(gs[:2,1]); mass=fig.add_subplot(gs[2,1])
    for i,(length,color) in enumerate(zip(lengths,colors)):
        x=np.arange(len(signed[length]))-(len(signed[length])-1)/2
        spans.plot([-(length-1)/2,(length-1)/2],[4-i,4-i],color=color,lw=3)
        spans.text(105,4-i,f"{length} bp",color=color,va='center',fontsize=10)
        wave.plot(x,signed[length],color=color,lw=2,ls='--' if length>mode else '-',label=f"{length} bp")
    spans.set_xlim(-110,110);spans.set_ylim(-0.7,4.7);spans.axis('off');spans.set_title('Observed fragments aligned at their centres',fontsize=11)
    wave.set(xlabel='Position from fragment centre (bp)',ylabel='PNS contribution per base',xlim=(-110,110))
    wave.axhline(0,color='0.4',lw=.8);wave.axvline(0,color='0.5',lw=.8,ls=':')
    wave.legend(ncol=3,fontsize=9,loc='lower center',frameon=False)
    lengths_all=np.arange(100,235)
    all_signed,_=precompute_distributions(lengths_all,mode)
    amplitude.plot(lengths_all,[all_signed[x].max() for x in lengths_all],color='#cf3e4e',lw=2)
    amplitude.axvline(mode,color='0.5',lw=.8,ls=':')
    amplitude.set(xlabel='Fragment length (bp)',ylabel='Maximum PNS contribution',title='Peak amplitude is highest at the mode')
    mass.plot(lengths_all,[all_signed[x][all_signed[x]>0].sum() for x in lengths_all],color='#cf3e4e',label='Positive mass')
    mass.plot(lengths_all,[all_signed[x][all_signed[x]<0].sum() for x in lengths_all],color='#377bb5',label='Negative mass')
    mass.set(xlabel='Fragment length (bp)',ylabel='Signed mass',ylim=(-135,135),yticks=[-100,0,100])
    mass.legend(fontsize=9,loc='center',ncol=2,frameon=False)
    for ax in [wave,amplitude,mass]:
        ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.2)
    fig.suptitle('PNS length adaptation · equal mass, changing width and amplitude',fontsize=14,weight='bold')
    save(fig,directory,'pns_length_adaptation_mode167')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parents[1]/'docs/images')
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','savefig.facecolor':'white'})
    geometry_figure(args.output_dir);adaptation_figure(args.output_dir)


if __name__=='__main__':
    main()
