"""
parse bonds from reaxff/bonds output for polyester/HMMA crosslinking simulation

crosslinking reaction:
    polyester-OH + HMMA-N-CH2-OCH3 -> polyester-O-CH2-N + CH3OH

three things are tracked:
    1. crosslinked oxygens: OH groups that have lost their H and gained a C bond
    2. remaining OH groups: unreacted hydroxyl sites
    3. methanol formed: CH3OH leaving group confirms transetherification occurred

bonds.reaxff column format:
    id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q

atom types (specorder C H O N):
    1 = C, 2 = H, 3 = O, 4 = N

usage:
    python parse_bonds_polyester.py
    python parse_bonds_polyester.py bonds.reaxff
    python parse_bonds_polyester.py bonds.reaxff --bo-cutoff 0.3
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
                ts      = int(line.split()[2])
                current = {'timestep': ts, 'atoms': {}}
                continue

            if line.startswith('#'):
                continue

            if current is None:
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


def find_hydroxyl_oxygens(snapshot, bo_cutoff=0.5):
    """
    find all oxygen atoms bonded to hydrogen with bond order above cutoff
    these are the reactive OH sites at the start of the simulation
    returns a set of oxygen atom ids
    """
    hydroxyl_oxygens = set()
    atoms = snapshot['atoms']

    for atom_id, data in atoms.items():
        if data['type'] != 3:   # skip non-oxygen atoms
            continue
        for partner_id, bo in data['partners']:
            partner_type = atoms.get(partner_id, {}).get('type', 0)
            if partner_type == 2 and bo >= bo_cutoff:   # oxygen bonded to hydrogen
                hydroxyl_oxygens.add(atom_id)
                break

    return hydroxyl_oxygens


def count_crosslinks(snapshots, bo_cutoff=0.3, oh_loss_cutoff=0.3):
    """
    at each timestep count:
        - how many original OH oxygens have crosslinked (lost H, gained C)
        - how many OH groups remain unreacted
        - how many methanol molecules have formed (CH3OH leaving group)

    a crosslink is defined as an oxygen that:
        1. had an O-H bond at t=0
        2. no longer has an O-H bond (bo < oh_loss_cutoff)
        3. now has a new C-O bond (bo > bo_cutoff)
    """

    if not snapshots:
        print("error: no snapshots found")
        sys.exit(1)

    # identify reactive OH sites from the first snapshot
    initial_oh   = find_hydroxyl_oxygens(snapshots[0], bo_cutoff=0.5)
    n_initial_oh = len(initial_oh)
    print(f"\n  initial hydroxyl oxygens (reactive sites): {n_initial_oh}")

    if n_initial_oh == 0:
        types  = [d['type'] for d in snapshots[0]['atoms'].values()]
        counts = Counter(types)
        print(f"  atom type distribution: {dict(counts)}")
        print(f"  expected: 1=C, 2=H, 3=O, 4=N")
        for atom_id, data in snapshots[0]['atoms'].items():
            if data['type'] == 3:
                print(f"  sample O atom {atom_id}: partners = {data['partners'][:3]}")
                break

    timesteps        = []
    crosslink_counts = []
    oh_remaining     = []
    methanol_counts  = []

    for snap in snapshots:
        ts    = snap['timestep']
        atoms = snap['atoms']

        crosslinked = 0
        still_oh    = 0
        methanol    = 0

        # check each original hydroxyl oxygen
        for o_id in initial_oh:
            if o_id not in atoms:
                continue

            has_h_bond = False
            has_c_bond = False

            for partner_id, bo in atoms[o_id]['partners']:
                ptype = atoms.get(partner_id, {}).get('type', 0)
                if ptype == 2 and bo >= oh_loss_cutoff:
                    has_h_bond = True
                if ptype == 1 and bo >= bo_cutoff:
                    has_c_bond = True

            if has_h_bond:
                still_oh += 1
            elif has_c_bond:
                crosslinked += 1

        # count methanol molecules: carbon with 0 C-C bonds, 1 C-O bond, 3 C-H bonds
        # where the bonded oxygen also has an O-H bond
        for atom_id, data in atoms.items():
            if data['type'] != 1:
                continue

            c_bonds = [(pid, bo) for pid, bo in data['partners']
                       if atoms.get(pid, {}).get('type', 0) == 1 and bo >= bo_cutoff]
            o_bonds = [(pid, bo) for pid, bo in data['partners']
                       if atoms.get(pid, {}).get('type', 0) == 3 and bo >= bo_cutoff]
            h_bonds = [(pid, bo) for pid, bo in data['partners']
                       if atoms.get(pid, {}).get('type', 0) == 2 and bo >= bo_cutoff]

            if len(c_bonds) == 0 and len(o_bonds) == 1 and len(h_bonds) == 3:
                o_id = o_bonds[0][0]
                if o_id in atoms:
                    o_h = [(pid, bo) for pid, bo in atoms[o_id]['partners']
                           if atoms.get(pid, {}).get('type', 0) == 2
                           and bo >= bo_cutoff]
                    if len(o_h) >= 1:
                        methanol += 1

        timesteps.append(ts)
        crosslink_counts.append(crosslinked)
        oh_remaining.append(still_oh)
        methanol_counts.append(methanol)

    return (np.array(timesteps), np.array(crosslink_counts),
            np.array(oh_remaining), np.array(methanol_counts),
            n_initial_oh)


def plot_results(timesteps, crosslink_counts, oh_remaining,
                 methanol_counts, n_initial_oh,
                 timestep_fs=0.25, output='polyester_crosslinks.png'):
    """plot crosslink conversion, remaining OH and methanol formation vs time"""

    time_ps    = timesteps * timestep_fs / 1000.0
    conversion = (crosslink_counts / n_initial_oh * 100.0
                  if n_initial_oh > 0
                  else np.zeros_like(crosslink_counts, dtype=float))

    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)

    ax1 = axes[0]
    ax1.plot(time_ps, crosslink_counts, 'b-', lw=1.5, label='crosslinked O', alpha=0.8)
    ax1.plot(time_ps, oh_remaining,     'r-', lw=1.5, label='remaining OH',  alpha=0.8)
    ax1.axhline(y=n_initial_oh, color='gray', ls='--', alpha=0.5,
                label=f'initial OH ({n_initial_oh})')
    ax1.set_ylabel('count')
    ax1.set_title('polyester/HMMA crosslinking — ReaxFF')
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
    ax3.plot(time_ps, methanol_counts, 'm-', lw=1.5, alpha=0.8)
    ax3.set_ylabel('methanol molecules')
    ax3.set_xlabel('simulation time (ps)')
    ax3.set_title('CH3OH leaving group — confirms transetherification')
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output, bbox_inches='tight', dpi=150)
    plt.close()
    print(f"\nsaved: {output}")


def print_summary(timesteps, crosslink_counts, oh_remaining,
                  methanol_counts, n_initial_oh, timestep_fs=0.25):
    """print summary statistics to terminal"""

    time_ps          = timesteps * timestep_fs / 1000.0
    final_crosslinks = int(crosslink_counts[-1]) if len(crosslink_counts) > 0 else 0
    final_oh         = int(oh_remaining[-1])     if len(oh_remaining)     > 0 else 0
    final_methanol   = int(methanol_counts[-1])  if len(methanol_counts)  > 0 else 0
    final_conv       = (final_crosslinks / n_initial_oh * 100) if n_initial_oh > 0 else 0
    total_time       = time_ps[-1]               if len(time_ps)          > 0 else 0

    print(f"\n  total simulation time:       {total_time:.1f} ps")
    print(f"  initial OH groups:           {n_initial_oh}")
    print(f"  final crosslinked oxygens:   {final_crosslinks}")
    print(f"  remaining OH groups:         {final_oh}")
    print(f"  crosslink conversion:        {final_conv:.1f}%")
    print(f"  methanol molecules formed:   {final_methanol}")
    print(f"  mass balance check:          {final_crosslinks + final_oh} / {n_initial_oh}")

    if len(crosslink_counts) > 10 and n_initial_oh > 0:
        half_max = n_initial_oh * 0.5
        idx_half = np.where(crosslink_counts >= half_max)[0]
        if len(idx_half) > 0:
            print(f"  time to 50% conversion:      {time_ps[idx_half[0]]:.1f} ps")
        else:
            print(f"  50% conversion not reached in simulation window")


# parse arguments and run
parser = argparse.ArgumentParser(
    description='analyse polyester/HMMA crosslinking from reaxff/bonds output')
parser.add_argument('bonds_file', nargs='?', default='bonds.reaxff',
                    help='reaxff/bonds output file (default: bonds.reaxff)')
parser.add_argument('--bo-cutoff', type=float, default=0.3)
parser.add_argument('--timestep',  type=float, default=0.25)
parser.add_argument('--output',    default='polyester_crosslinks.png')
args = parser.parse_args()

if not os.path.exists(args.bonds_file):
    print(f"error: cannot find {args.bonds_file}")
    sys.exit(1)

snapshots = parse_bonds_file(args.bonds_file)

(timesteps, crosslink_counts,
 oh_remaining, methanol_counts,
 n_initial_oh) = count_crosslinks(snapshots, bo_cutoff=args.bo_cutoff)

print_summary(timesteps, crosslink_counts, oh_remaining,
              methanol_counts, n_initial_oh, timestep_fs=args.timestep)

if n_initial_oh > 0:
    plot_results(timesteps, crosslink_counts, oh_remaining,
                 methanol_counts, n_initial_oh,
                 timestep_fs=args.timestep, output=args.output)
else:
    print("\nno OH groups found -- check atom type mapping")
    print("run: grep 'Masses' system_polyester.data -A 6")
