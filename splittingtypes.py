"""
This script produces the table of splitting types and their associated
decomposition and inertia groups.

Will eventually work generally for any group, but for now, only works
works for the wild groups over 2.

@author Kevin H. Wilson
"""
from __future__ import print_function


# Sage for wild case
G = DihedralGroup(4)
sigma = G((1, 2, 3, 4))
tau = G([(1, 4),(2, 3)])
V41 = G.subgroup([tau, sigma^2])
V42 = G.subgroup([sigma * tau, sigma^2])
D4 = G

for D_name, D in [('V41', V41), ('V42', V42), ('D4', D4)]:
    print("---------- D = {} ----------".format(D_name))
    for H_name, H in [('M', G.subgroup([])),
                      ('L1', G.subgroup([tau])), ('K1', G.subgroup([tau, sigma^2])),
                      ('L2', G.subgroup([sigma * tau])), ('K2', G.subgroup([sigma * tau, sigma^2])),
                      ('L', G.subgroup([sigma^2])), ('K', G.subgroup([sigma]))]:
        cosets = G.cosets(H)
        orbits = set()
        for coset in cosets:
            orbits.add(frozenset(d * c for d in D for c in coset))

        print(H_name, [len(orbit) / len(H) for orbit in orbits])
    print()
