import click
from sympy.combinatorics import Permutation

from . import expected_number, splitting_types
from .splitting_types import GaloisGroup


@click.group()
def cli():
    """ Computations for counting D4 fields """


@cli.command("expectation")
def expectation_command():
    """
    Compute Table 2: the expectated number of fields at p = 2
    """
    outcomes = expected_number.compute()
    assert len(outcomes) == 14
    print("d = 0")
    for q in [0, 2, 3, 4, 5, 6]:
        print("q = {}: {}".format(q, outcomes[(0, q)]))
    print("d = 2")
    for q in [0, 2, 4, 5]:
        print("q = {}: {}".format(q, outcomes[(2, q)]))
    print("d = 3")
    for q in [0, 2, 4, 5]:
        print("q = {}: {}".format(q, outcomes[(3, q)]))
    print(
        "Total weight:", sum(val / 2 ** (d + q + 3) for (d, q), val in outcomes.items())
    )
    print()
    print("Latex table:")
    expected_number.print_table(outcomes)


@cli.command("splitting")
def splitting_command():
    """
    Compute Table 1: splitting types associated with D4 fields
    """
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
    # Conjugate of Gal_M_L1, so unused in computation
    # Gal_M_L1_prime = GaloisGroup(
    #     group=(ι, σ ** 2 * τ), latex_name=r"\langle \sigma^2 \tau \rangle"
    # )

    Gal_M_L2 = GaloisGroup(group=(ι, σ * τ), latex_name=r"\langle \sigma \tau \rangle")
    # Conjugate of Gal_M_L2, so unused in computation
    # Gal_M_L2_prime = GaloisGroup(
    #     group=(ι, σ ** 3 * τ), latex_name=r"\langle \sigma^3 \tau \rangle"
    # )

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
        Gal_M_K1,
        Gal_M_L2,
        Gal_M_K2,
        Gal_M_L,
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

    # Setup the header
    print(r"\begin{tabular}{|c|c|||c||c|c||c|c||c|c|||c|c|}")
    print(r"  \hline &&&&&&&&&\multicolumn{2}{|c|}{}\\")

    field_names = ["M", "L_1", "K_1", "L_2", "K_2", "L", "K"]
    print(
        r"  $I_p$ & $D_p$ & "
        + " & ".join(f"$\\varsigma_p({name})$" for name in field_names)
        + r" & \multicolumn{2}{|c|}{} \\[12pt]"
    )
    print(r"  \hline &&&&&&&&&&\\[-10.5pt]")
    print(r"  \hline &&&&&&&&&&\\[-10.5pt]")
    print(r"  \hline &&&&&&&&&&\\[-5pt]")
    print()

    # Begin writing out the details
    for row_num, (inertia_group, decomposition_group) in enumerate(
        zip(inertia_groups, decomposition_groups)
    ):
        print("  ", end="")
        print(f"${inertia_group.latex_name}$ & ${decomposition_group.latex_name}$ & ")
        to_print = []
        for galois_group in gal_groups:
            cosets = splitting_types.get_cosets(D4, galois_group)
            orbits = splitting_types.get_orbits(cosets, decomposition_group)
            inertia = splitting_types.compute_inertia(orbits, inertia_group)
            to_print.append(f"${splitting_types.write_inertia(inertia)}$")

        print("  " + " & ".join(to_print), end="")

        # Deal with the niceties of the labels
        if row_num == 0:
            print(
                r" & \multirow{5}{*}{\rotatebox[origin=c]{270}{\small Unramified}} & \multirow{10}{*}{\rotatebox[origin=c]{270}{\small Lacks central inertia}}\\"
            )
        elif row_num == 4:
            print(r" & & \\")
            print(r"    [5pt] \cline{1-10} &&&&&&&&&&\\[-4pt]")
        elif row_num == 5:
            print(r" & \multirow{4}{*}{\rotatebox[origin=c]{270}{\small Tame}} & \\")
        elif row_num == 8:
            print(r" & & \\")
            print(r"    [5pt] \hline &&&&&&&&&&\\[-10.5pt]\hline &&&&&&&&&&\\[-4pt]")
        elif row_num == 9:
            print(
                r" & \multirow{6}{*}{\rotatebox[origin=c]{270}{\small Tame}} & \multirow{10}{*}{\rotatebox[origin=c]{270}{\small Has central inertia}}\\"
            )
        elif row_num == 14:
            print(r" & & \\")
            print(r"    [5pt] \cline{1-10} &&&&&&&&&&\\[-4pt]")
        elif row_num == 15:
            print(r" & \multirow{3}{*}{\rotatebox[origin=c]{270}{\small Wild}} & \\")
        else:
            print(r" & & \\")

        print()

    # Close out the table
    print(r"  \hline")
    print(r"\end{tabular}\\[.1in]")
