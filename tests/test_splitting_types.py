from sympy.combinatorics import Permutation

from d4counting import splitting_types
from d4counting.splitting_types import GaloisGroup


def test_get_cosets():
    """
    Try a few simple examples with cyclic groups
    """
    σ = Permutation(0, 1, 2, 3)
    ι = σ ** 4
    trivial_group = GaloisGroup(group=[ι], latex_name=r"\{ 1 \}")
    cyclic_2 = GaloisGroup(group=[ι, σ ** 2], latex_name=r"C_2")
    cyclic_4 = GaloisGroup(group=[ι, σ, σ ** 2, σ ** 3], latex_name=r"C_4")

    # First try the trival group; should be each element on its own
    cosets = splitting_types.get_cosets(cyclic_4, trivial_group)
    assert cosets == frozenset(frozenset({g}) for g in cyclic_4.group)

    # Next try C2 ⊆ C4. Should be two sets: C2 and the rest
    cosets = splitting_types.get_cosets(cyclic_4, cyclic_2)
    assert cosets == frozenset({frozenset(cyclic_2), frozenset({σ, σ ** 3})})

    # Finally, try the whole C4 ⊆ C$. Should be one set
    cosets = splitting_types.get_cosets(cyclic_4, cyclic_4)
    assert cosets == frozenset({frozenset(cyclic_4)})


def test_get_orbits():
    """
    Try a simple example with the cyclic group of order 4
    """
    σ = Permutation(0, 1, 2, 3)
    ι = σ ** 4
    trivial_group = GaloisGroup(group=[ι], latex_name=r"\{ 1 \}")
    cyclic_2 = GaloisGroup(group=[ι, σ ** 2], latex_name=r"C_2")
    cyclic_4 = GaloisGroup(group=[ι, σ, σ ** 2, σ ** 3], latex_name=r"C_4")

    # Look at trivial \\ C4 // C2
    cosets = splitting_types.get_cosets(cyclic_4, trivial_group)
    orbits = splitting_types.get_orbits(cosets, cyclic_2)
    assert orbits == frozenset(
        {
            frozenset({frozenset({ι}), frozenset({σ ** 2})}),
            frozenset({frozenset({σ}), frozenset({σ ** 3})}),
        }
    )


def test_compute_inertia():
    """
    Run a complete example for D4 that we know the answer to
    """
    σ = Permutation(0, 1, 2, 3)
    τ = Permutation(0, 1)(2, 3)
    ι = Permutation(3)  # The identity

    # Big Galois group
    Gal_M_Q = GaloisGroup(
        group=(ι, σ, σ ** 2, σ ** 3, τ, σ * τ, σ ** 2 * τ, σ ** 3 * τ),
        latex_name=r"D_4",
    )

    # Groups that fix the subfield
    Gal_M_L1 = GaloisGroup(group=(ι, τ), latex_name=r"\langle \tau \rangle")
    Gal_M_L2 = GaloisGroup(group=(ι, σ * τ), latex_name=r"\langle \tau \rangle")

    # The decomposition group
    Gal_M_K = GaloisGroup(
        group=(ι, σ, σ ** 2, σ ** 3), latex_name=r"\langle \sigma \rangle"
    )

    # The inertia group
    Gal_M_L = GaloisGroup(group=(ι, σ ** 2), latex_name=r"\langle \sigma^2 \rangle")

    cosets = splitting_types.get_cosets(Gal_M_Q, Gal_M_L1)
    orbits = splitting_types.get_orbits(cosets, Gal_M_K)
    inertia = splitting_types.compute_inertia(orbits, Gal_M_L)
    assert splitting_types.write_inertia(inertia) == "(2^2)"

    cosets = splitting_types.get_cosets(Gal_M_Q, Gal_M_L2)
    orbits = splitting_types.get_orbits(cosets, Gal_M_K)
    inertia = splitting_types.compute_inertia(orbits, Gal_M_L)
    assert splitting_types.write_inertia(inertia) == "(2^2)"

    cosets = splitting_types.get_cosets(Gal_M_Q, Gal_M_L)
    orbits = splitting_types.get_orbits(cosets, Gal_M_K)
    inertia = splitting_types.compute_inertia(orbits, Gal_M_L)
    assert splitting_types.write_inertia(inertia) == "(2 2)"
