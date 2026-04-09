"""
analysis script for DGEBF/DETDA epoxy crosslinking simulation

crosslinking reaction:
    DGEBF-epoxide + DETDA-NH2 -> DGEBF-CH(OH)-CH2-N-DETDA

three things are tracked:
    1. c-n bond formation: epoxide carbon bonding to amine nitrogen
    2. n-h bond loss: amine reactive sites consumed as reaction proceeds
    3. o-h bond formation: oxygen picking up proton after epoxide ring opening

atom types (specorder c h o n):
    1 = c, 2 = h, 3 = o, 4 = n

system: 50 dgebf x 2 epoxide groups = 100 potential crosslinks
        25 detda x 4 n-h bonds      = 100 potential crosslinks

usage:
    python epoxy_analysis.py
    python epoxy_analysis.py bonds.reaxff
    python epoxy_analysis.py bonds_hold.reaxff --output epoxy_hold.png
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
from collections import Counter


def parse_bonds_file(filename):
    """
    read the reaxff/bonds file and return a list of snapshots
    each snapshot contains every atom and its bond partners + bond orders at that timestep

    bonds.reaxff column format:
        id type nb id_1..id_nb mol bo_1..bo_nb abo nlp q
    """
    snapshots = []
    current   = None

    print(f"parsing {filename}...")

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # each snapshot starts with a timestep header
            if line.startswith('# Timestep'):
                if current is not None:
                    snapshots.append(current)
                current = {'timestep': int(line.split()[2]), 'atoms': {}}
                continue

            if line.startswith('#') or current is None:
                continue

            parts = line.split()
            if len(parts) < 4:
                continue

            try:
                atom_id   = int(parts[0])
                atom_type = int(parts[1])
                n_bonds   = int(parts[2])

                partners = []
                if n_bonds > 0:
                    # extract bond partner ids
                    bond_ids    = [int(parts[3 + i]) for i in range(n_bonds)]
                    # skip the mol column at position 3 + n_bonds
                    bo_start    = 3 + n_bonds + 1
                    # extract corresponding bond orders
                    bond_orders = [float(parts[bo_start + i]) for i in range(n_bonds)]
                    partners    = list(zip(bond_ids, bond_orders))

                current['atoms'][atom_id] = {
                    'type':     atom_type,
                    'partners': partners
                }

            except (ValueError, IndexError):
                continue

    if current is not None:
        snapshots.append(current)

    print(f"  loaded {len(snapshots)} snapshots")
    return snapshots


def count_crosslinks(snapshots, bo_cutoff=0.3, n_epoxide=100):
    """
    at each timestep count:
        - c-n bonds: crosslinks formed (counted from carbon side to avoid double counting)
        - n-h bonds: remaining amine reactive sites
        - o-h bonds: formed after epoxide ring opening, confirms correct mechanism

    n_epoxide: total epoxide groups in system, used to calculate conversion %
    """

    if not snapshots:
        print("error: no snapshots found")
        sys.exit(1)

    # diagnostic on first snapshot
    types  = [d['type'] for d in snapshots[0]['atoms'].values()]
    counts = Counter(types)
    print(f"\n  atom type distribution: {dict(counts)}")
    print(f"  (1=c, 2=h, 3=o, 4=n)")

    # count initial n-h bonds — these are the reactive sites
    initial_nh = 0
    for atom_id, data in snapshots[0]['atoms'].items():
        if data['type'] != 4:
            continue
        for partner_id, bo in data['partners']:
            ptype = snapshots[0]['atoms'].get(partner_id, {}).get('type', 0)
            if ptype == 2 and bo >= bo_cutoff:
                initial_nh += 1

    print(f"  initial n-h bonds (reactive sites): {initial_nh}")
    print(f"  theoretical max crosslinks (epoxide): {n_epoxide}")

    timesteps  = []
    cn_counts  = []
    nh_counts  = []
    oh_counts  = []

    for snap in snapshots:
        atoms = snap['atoms']

        cn = 0
        nh = 0
        oh = 0

        for atom_id, data in atoms.items():

            # count c-n bonds from carbon side only to avoid double counting
            if data['type'] == 1:
                for partner_id, bo in data['partners']:
                    if bo >= bo_cutoff:
                        if atoms.get(partner_id, {}).get('type', 0) == 4:
                            cn += 1

            # count remaining n-h bonds
            elif data['type'] == 4:
                for partner_id, bo in data['partners']:
                    if bo >= bo_cutoff:
                        if atoms.get(partner_id, {}).get('type', 0) == 2:
                            nh += 1

            # count o-h bonds formed after ring opening
            elif data['type'] == 3:
                for partner_id, bo in data['partners']:
                    if bo >= bo_cutoff:
                        if atoms.get(partner_id, {}).get('type', 0) == 2:
                            oh += 1

        timesteps.append(snap['timestep'])
        cn_counts.append(cn)
        nh_counts.append(nh)
        oh_counts.append(oh)

    return (np.array(timesteps), np.array(cn_counts),
            np.array(nh_counts), np.array(oh_counts),
            initial_nh, n_epoxide)


def plot_results(timesteps, cn_counts, nh_counts, oh_counts,
                 initial_nh, n_epoxide,
                 timestep_fs=0.25, output='epoxy_crosslinks.png'):
    """plot c-n bond count, conversion %, and mechanistic bond tracking vs time"""

    time_ps    = timesteps * timestep_fs / 1000.0
    conversion = cn_counts / n_epoxide * 100.0

    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)

    ax1 = axes[0]
    ax1.plot(time_ps, cn_counts, 'b-', lw=1.5, label='c-n bonds (crosslinks)', alpha=0.8)
    ax1.axhline(y=n_epoxide, color='r', ls='--', alpha=0.5,
                label=f'max crosslinks ({n_epoxide})')
    ax1.set_ylabel('c-n bond count')
    ax1.set_title('DGEBF/DETDA epoxy crosslinking — ReaxFF (Vashisth 2018 FF)')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    ax2 = axes[1]
    ax2.plot(time_ps, conversion, 'g-', lw=1.5, alpha=0.8)
    ax2.set_ylabel('crosslink conversion (%)')
    ax2.set_ylim(0, 105)
    ax2.grid(True, alpha=0.3)
    if len(conversion) > 0:
        fc = conversion[-1]
        ax2.annotate(f'final: {fc:.1f}%',
                     xy=(time_ps[-1], fc),
                     xytext=(time_ps[-1]*0.6, min(fc+10, 95)),
                     arrowprops=dict(arrowstyle='->', color='black'),
                     fontsize=11)

    ax3 = axes[2]
    ax3.plot(time_ps, nh_counts, 'r-',  lw=1.5, label='n-h bonds remaining', alpha=0.8)
    ax3.plot(time_ps, oh_counts, 'c--', lw=1.5, label='o-h bonds (ring opening)', alpha=0.8)
    ax3.axhline(y=initial_nh, color='gray', ls='--', alpha=0.4,
                label=f'initial n-h ({initial_nh})')
    ax3.set_ylabel('bond count')
    ax3.set_xlabel('simulation time (ps)')
    ax3.set_title('reaction mechanism: n-h consumed, o-h formed')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output, bbox_inches='tight', dpi=150)
    plt.close()
    print(f"\nsaved: {output}")


def print_summary(timesteps, cn_counts, nh_counts, oh_counts,
                  initial_nh, n_epoxide, timestep_fs=0.25):
    """print summary statistics to terminal"""

    time_ps     = timesteps * timestep_fs / 1000.0
    final_cn    = int(cn_counts[-1])  if len(cn_counts)  > 0 else 0
    final_nh    = int(nh_counts[-1])  if len(nh_counts)  > 0 else 0
    final_oh    = int(oh_counts[-1])  if len(oh_counts)  > 0 else 0
    final_conv  = final_cn / n_epoxide * 100 if n_epoxide > 0 else 0
    total_time  = time_ps[-1]         if len(time_ps)    > 0 else 0
    nh_consumed = initial_nh - final_nh

    print(f"\n  total simulation time:       {total_time:.1f} ps")
    print(f"  epoxide groups (max xl):     {n_epoxide}")
    print(f"  initial n-h bonds:           {initial_nh}")
    print(f"  final c-n bonds (crosslinks):{final_cn}")
    print(f"  crosslink conversion:        {final_conv:.1f}%")
    print(f"  n-h bonds consumed:          {nh_consumed}")
    print(f"  o-h bonds formed:            {final_oh}")
    print(f"  mechanistic check (cn=oh?):  cn={final_cn}, oh change={final_oh}")
    print(f"  vashisth 2018 experimental:  ~85-95%")

    if len(cn_counts) > 10:
        half_max = n_epoxide * 0.5
        idx_half = np.where(cn_counts >= half_max)[0]
        if len(idx_half) > 0:
            print(f"  time to 50% conversion:      {time_ps[idx_half[0]]:.1f} ps")
        else:
            print(f"  50% conversion not reached in simulation window")


# parse arguments and run
parser = argparse.ArgumentParser(
    description='analyse DGEBF/DETDA crosslinking from reaxff/bonds output')
parser.add_argument('bonds_file', nargs='?', default='bonds.reaxff',
                    help='reaxff/bonds file (default: bonds.reaxff)')
parser.add_argument('--bo-cutoff',  type=float, default=0.3)
parser.add_argument('--n-epoxide',  type=int,   default=100)
parser.add_argument('--timestep',   type=float, default=0.25)
parser.add_argument('--output',     default='epoxy_crosslinks.png')
args = parser.parse_args()

if not os.path.exists(args.bonds_file):
    print(f"error: cannot find {args.bonds_file}")
    sys.exit(1)

snapshots = parse_bonds_file(args.bonds_file)

(timesteps, cn_counts, nh_counts,
 oh_counts, initial_nh,
 n_epoxide) = count_crosslinks(snapshots,
                               bo_cutoff=args.bo_cutoff,
                               n_epoxide=args.n_epoxide)

print_summary(timesteps, cn_counts, nh_counts, oh_counts,
              initial_nh, n_epoxide, timestep_fs=args.timestep)

plot_results(timesteps, cn_counts, nh_counts, oh_counts,
             initial_nh, n_epoxide,
             timestep_fs=args.timestep, output=args.output)
