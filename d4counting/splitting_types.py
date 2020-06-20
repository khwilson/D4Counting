from dataclasses import dataclass
from typing import List, FrozenSet, Tuple

from sympy.combinatorics import Permutation


SetOfCosets = FrozenSet[FrozenSet[Permutation]]


@dataclass(eq=True, frozen=True)
class GaloisGroup:
    """ A simple representation of a group that also includes a name for printing purposes """

    group: Tuple[Permutation]
    latex_name: str

    def __iter__(self):
        return iter(self.group)


def get_cosets(big_galois: GaloisGroup, small_galois: GaloisGroup) -> SetOfCosets:
    """
    Given a big group `big_galois` and a subgroup `small_galois`, return the cosets
    of small_galois \\ big_galois.

    Args:
        big_galois: A `GaloisGroup` whose cosets to examine
        small_galois: The acting subgroup of big_galois

    Returns:
        A colletion of cosets. Each coset is a frozenset of `Permutation`s.
    """
    return frozenset(frozenset(h * g for h in small_galois) for g in big_galois)


def get_orbits(cosets: SetOfCosets, group: GaloisGroup) -> FrozenSet[SetOfCosets]:
    """
    Given a collection of cosets, ask for the orbits of `group` when acting on those
    cosets. Note that there are no checks here to make sure, e.g., the action of
    `group` makes sense or that `cosets` is closed under this action. Caveat Emptor.

    Args:
        cosets: The cosets to examine
        group: The group whose action to compute

    Returns:
        A collection of collections of cosets. Each collection represents an orbit.
        Each coset is represented as frozenset.
    """
    return frozenset(
        frozenset(frozenset(h * g for h in coset) for g in group) for coset in cosets
    )


def compute_inertia(
    decomp_orbits: FrozenSet[SetOfCosets], inertia_group: GaloisGroup
) -> List[Tuple[int, int]]:
    """
    Given a colletion of orbits of the decomposition group (as returned by `get_orbits`),
    compute the inertial decomposition for the given inertia group.

    Args:
        decomp_orbits: The orbits of the decomposition group on cosets of small_galois \\ big_galois
        inertia_group: The inertia group to examine

    Returns:
        The inertial decomposition. Each element is a pair (f, e) where:
            * f is the inertia degree of the associated prime and
            * e is the ramification index
    """
    inertia = []
    for orbit in decomp_orbits:
        suborbits = get_orbits(orbit, inertia_group)
        example_orbit = next(iter(suborbits))
        inertia.append((len(suborbits), len(example_orbit)))
    return inertia


def write_inertia(inertia: List[Tuple[int, int]]):
    """
    Given an inertia decomposition as returned by `compute_inertia`, create the string
    version that is latex displayable.

    Example:
        >>> write_inertia([(2, 1), (2, 3), (1, 2), (1, 1)])
        '(1^2 1 2^3 2)'

    Args:
        inertia: The inertia representation as returned by `compute_inertia`
    """
    inertia.sort(key=lambda x: (x[0], -x[1]))
    elts = [f"{f}^{e}" if e > 1 else f"{f}" for f, e in inertia]
    return "(" + " ".join(elts) + ")"


def main():
    σ = Permutation(0, 1, 2, 3)
    τ = Permutation(0, 1)(2, 3)
    ι = Permutation(3)  # The identity

    # Make sure these are actually generators of D4
    assert τ ** 2 == ι
    assert σ ** 4 == ι
    assert τ * σ * τ == σ ** 3

    # The big group
    D4 = GaloisGroup(
        group=(ι, σ, σ ** 2, σ ** 3, τ, σ * τ, σ ** 2 * τ, σ ** 3 * τ),
        latex_name=r"D_4",
    )

    # The individual Galois groups corresponding to the field diagram in the paper
    Gal_M_M = GaloisGroup(group=(ι,), latex_name=r"\{1\}")
    Gal_M_L1 = GaloisGroup(group=(ι, τ), latex_name=r"\langle \tau \rangle")
    Gal_M_L1_prime = GaloisGroup(
        group=(ι, σ ** 2 * τ), latex_name=r"\langle \sigma^2 \tau \rangle"
    )
    Gal_M_L2 = GaloisGroup(group=(ι, σ * τ), latex_name=r"\langle \sigma \tau \rangle")
    Gal_M_L2_prime = GaloisGroup(
        group=(ι, σ ** 3 * τ), latex_name=r"\langle \sigma^3 \tau \rangle"
    )
    Gal_M_L = GaloisGroup(group=(ι, σ ** 2), latex_name=r"\langle \sigma^2 \rangle")
    Gal_M_K1 = GaloisGroup(
        group=(ι, τ, σ ** 2, σ ** 2 * τ), latex_name=r"\langle \tau, \sigma^2 \rangle"
    )
    Gal_M_K2 = GaloisGroup(
        group=(ι, σ * τ, σ ** 2, σ ** 3 * τ),
        latex_name=r"\langle \sigma \tau, \sigma^2 \rangle",
    )
    Gal_M_K = GaloisGroup(
        group=(ι, σ, σ ** 2, σ ** 3), latex_name=r"\langle \sigma \rangle"
    )

    gal_groups = [
        Gal_M_M,
        Gal_M_L1,
        Gal_M_L1_prime,
        Gal_M_L2,
        Gal_M_L2_prime,
        Gal_M_L,
        Gal_M_K1,
        Gal_M_K2,
        Gal_M_K,
    ]

    inertia_groups = [
        # Unramified
        Gal_M_M,
        Gal_M_M,
        Gal_M_M,
        Gal_M_M,
        Gal_M_M,
        # Tame, no central inertia
        Gal_M_L1,
        Gal_M_L1,
        Gal_M_L2,
        Gal_M_L2,
        # Tame, Central Inertia
        Gal_M_K,
        Gal_M_K,
        Gal_M_L,
        Gal_M_L,
        Gal_M_L,
        Gal_M_L,
        # Wild, Central Inertia
        Gal_M_K1,
        Gal_M_K2,
        D4,
    ]

    decomposition_groups = [
        # Unramified
        Gal_M_M,
        Gal_M_L,
        Gal_M_L2,
        Gal_M_L1,
        Gal_M_K,
        # Tame, no central inertia
        Gal_M_L1,
        Gal_M_K1,
        Gal_M_L2,
        Gal_M_K2,
        # Tame, Central Inertia
        Gal_M_K,
        D4,
        Gal_M_L,
        Gal_M_K1,
        Gal_M_K2,
        Gal_M_K,
        # Wild, Central Inertia
        Gal_M_K1,
        Gal_M_K2,
        D4,
    ]
