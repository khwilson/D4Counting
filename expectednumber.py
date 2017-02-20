"""
This file is used to compute the expected number of fields with
invariants (q, d) over 2

@author Kevin H. Wilson
"""
from __future__ import division, print_function

from collections import defaultdict, namedtuple
from io import StringIO
import csv
import sys
import textwrap


def get_quad_disc(disc):
    if disc == '*':
        return 0
    if '2' in disc:
        return 3
    assert '-' in disc
    return 2


C2Field = namedtuple('C2Field', 'c e f d eps poly G inertia slopes')
c2_fields_raw = textwrap.dedent("""\
    c,e,f,d,eps,Polynomial,G,inertia,slopes
    0,1,2,*,1,x2-x+1,C2,unram,[]
    2,2,1,-1,i,x2-2x+2,C2,C2,[2]
    2,2,1,-*,i,x2-2x+6,C2,C2,[2]
    3,2,1,-2,i,x2+2,C2,C2,[3]
    3,2,1,2*,-1,x2+6,C2,C2,[3]
    3,2,1,2,i,x2+10,C2,C2,[3]
    3,2,1,-2*,-1,x2+14,C2,C2,[3]
""")


C4Field = namedtuple('C4Field', 'c e f d eps poly G inertia slopes deg2_subfield')
c4_fields_raw=textwrap.dedent("""\
    c,e,f,d,eps,poly,G,inertia,slopes,deg2_subfield
    0,1,4,*,1,x4-x+1,C4,unram,"[]",*
    4,2,2,*,-1,x4-x2+5,C4,C2,"[2]",*
    6,2,2,*,1,x4+2x2+20,C4,C2,"[3]",*
    6,2,2,*,-1,x4-2x2+20,C4,C2,"[3]",*
    11,4,1,2,1,x4+12x2+2,C4,C4,"[3, 4]",2
    11,4,1,2,-1,x4+8x+14,C4,C4,"[3, 4]",2
    11,4,1,2,-1,x4+4x2+18,C4,C4,"[3, 4]",2
    11,4,1,2,1,x4+12x2+18,C4,C4,"[3, 4]",2
    11,4,1,2*,1,x4+8x2+8x+22,C4,C4,"[3, 4]",2*
    11,4,1,2*,-1,x4+4x2+10,C4,C4,"[3, 4]",2*
    11,4,1,2*,1,x4+12x2+10,C4,C4,"[3, 4]",2*
    11,4,1,2*,-1,x4+8x+6,C4,C4,"[3, 4]",2*
""")


V4Field = namedtuple('V4Field', 'c e f d eps poly G inertia slopes')
v4_fields_raw = textwrap.dedent("""\
    4,2,2,1,-1,x4+8x2+4,V4,C2,"[2]"
    6,2,2,1,-1,x4-6x2+4,V4,C2,"[3]"
    6,2,2,1,1,x4-2x2+4,V4,C2,"[3]"
    8,4,1,1,-1,x4+2x2+4x+10,V4,V4,"[2, 3]"
    8,4,1,1,-1,x4+6x2+1,V4,V4,"[2, 3]"
    8,4,1,1,1,x4+6x2+4x+14,V4,V4,"[2, 3]"
    8,4,1,1,1,x4+6x2+4x+6,V4,V4,"[2, 3]"
""")


D4Field = namedtuple('D4Field', 'c e f d eps poly G inertia slopes deg2_subfield')
d4_fields_raw = textwrap.dedent("""\
    c,e,f,d,eps,Polynomial,G,inertia,slopes,deg2_subfield
    4,2,2,-1,-i,x4+2x2+4x+4,D4,V4,"[2, 2]",*
    4,2,2,-*,-i,x4-5,D4,V4,"[2, 2]",*
    6,2,2,-1,i,x4+2x2-4,D4,V4,"[2,3]",*
    6,2,2,-*,-i,x4-20,D4,V4,"[2, 3]",*
    6,4,1,*,1,x4+2x3+2,D4,V4,"[2, 2]",-1
    6,4,1,*,1,x4+2x3+6,D4,V4,"[2, 2]",-*
    8,4,1,*,1,x4+2x2+4x+6,D4,V4,"[2, 3]",-*
    8,4,1,*,-1,x4+6x2+4x+2,D4,V4,"[2, 3]",-1
    9,4,1,2,1,x4+6x2+2,D4,D4,"[2, 3, 7/2]",-1
    9,4,1,2,-1,x4-2x2+2,D4,D4,"[2, 3, 7/2]",-1
    9,4,1,2*,1,x4+6x2+10,D4,D4,"[2, 3, 7/2]",-1
    9,4,1,2*,-1,x4+2x2+10,D4,D4,"[2, 3, 7/2]",-1
    9,4,1,-2,i,x4+2x2-2,D4,D4,"[2, 3, 7/2]",-*
    9,4,1,-2,-i,x4-2x2-2,D4,D4,"[2, 3, 7/2]",-*
    9,4,1,-2*,i,x4+2x2+6,D4,D4,"[2, 3, 7/2]",-*
    9,4,1,-2*,-i,x4-2x2+6,D4,D4,"[2, 3, 7/2]",-*
    10,4,1,-1,i,x4+2x2-9,D4,D4,"[2, 3, 7/2]",2*
    10,4,1,-1,i,x4+2x2-1,D4,D4,"[2, 3, 7/2]",2
    10,4,1,-1,-i,x4+6x2-9,D4,D4,"[2, 3, 7/2]",2
    10,4,1,-1,-i,x4+6x2-1,D4,D4,"[2, 3, 7/2]",2*
    10,4,1,-*,i,x4-6x2+3,D4,D4,"[2, 3, 7/2]",-2*
    10,4,1,-*,-i,x4+6x2+3,D4,D4,"[2, 3, 7/2]",-2*
    10,4,1,-*,i,x4-2x2+3,D4,D4,"[2, 3, 7/2]",-2
    10,4,1,-*,-i,x4+2x2+3,D4,D4,"[2, 3, 7/2]",-2
    11,4,1,2,-1,x4+2,D4,D4,"[2, 3, 4]",-2
    11,4,1,2,-1,x4+18,D4,D4,"[2, 3, 4]",-2
    11,4,1,2*,-1,x4+10,D4,D4,"[2, 3, 4]",-2*
    11,4,1,2*,-1,x4+26,D4,D4,"[2, 3, 4]",-2*
    11,4,1,-2,-i,x4+4x2+14,D4,C4,"[3, 4]",-2*
    11,4,1,-2,i,x4+8x+10,D4,C4,"[3, 4]",2*
    11,4,1,-2,i,x4+30,D4,D4,"[2, 3, 4]",2
    11,4,1,-2,i,x4+14,D4,D4,"[2, 3, 4]",2
    11,4,1,-2*,-i,x4+4x2+6,D4,C4,"[3, 4]",-2
    11,4,1,-2*,i,x4+12x2+6,D4,C4,"[3, 4]",-2
    11,4,1,-2*,i,x4+22,D4,D4,"[2, 3, 4]",2*
    11,4,1,-2*,i,x4+6,D4,D4,"[2, 3, 4]",2*
""")


def make_c2_fields():
    return parse_raw(c2_fields_raw, C2Field)

def make_c4_fields():
    return parse_raw(c4_fields_raw, C4Field)

def make_v4_fields():
    return parse_raw(v4_fields_raw, V4Field)

def make_d4_fields():
    return parse_raw(d4_fields_raw, D4Field)


def do_unram_fields():
    return [(0, 0, 8)]


def do_c2_fields():
    output = []
    for field in make_c2_fields():
        if not field.inertia == 'unram':
            output.extend([(0, field.c, 2),  # I = tau
                           (field.c, field.c, 2),  # I = sigma*tau
                           (0, 2 * field.c, 1) # I = sigma^2
                          ])
    return output

def do_c4_fields():
    output = []
#    return output
    for field in make_c4_fields():
        if field.inertia == 'unram':
            continue
        elif field.inertia == 'C2':
            output.append((0, field.c, 2))  # 2 for automorphisms
        else:
            # I = sigma so survives dividing by <tau, sigma^2>
            # Thus, the quadratic subfield of the field is d
            d = get_quad_disc(field.deg2_subfield)

            # chi(<sigma>) = 0 and chi(<sigma^2>) = 0 so the conductor is
            # 2 * slope[0] + 2 * (slope[1] - slope[0]) == 2 * slope[1]
            c = 2 * field.slopes[-1]

            # And the number of fields is just the number of automorphisms
            output.append((d, c, 2))
    return output


def do_v4_fields():
    output = []
    for field in make_v4_fields():
        if field.inertia == 'unram':
            raise Exception("Not possible")

        elif field.inertia == 'C2':
            # If V4 = <tau, sigma^2>, then inertia never survives and so d = 0
            # On the other hand, the conductor depends on whether the image of
            # inertia is <sigma^2> or a reflection. If it is a reflection, then
            # the conductor is just the slope. If it is the center, then
            # it is *twice* the slope
            output.extend([(0, field.slopes[-1], 4),  # 4/6 automorphisms yield reflections
                           (0, 2 * field.slopes[-1], 2)
                          ])

            # If V4 = <sigma*tau, sigma^2> then inertia survives when it is a reflection,
            # but otherwise dies
            output.extend([(field.slopes[-1], field.slopes[-1], 4),  # 4/6 automorphisms yield reflections
                           (0, 2 * field.slopes[-1], 2)
                          ])

        elif field.inertia == 'V4':
#            continue
            # If V4 = <tau, sigma^2>, then inertia never survives so d = 0
            # On the other hand, the conductor depends on the image of
            # the *second* inertia group.
            output.extend([(0, 2 * field.slopes[0] + (field.slopes[1] - field.slopes[0]), 4),
                           (0, 2 * field.slopes[1], 2)
                          ])

            # If V4 = <sigma*tau, sigma^2>, then inertia always survives. However,
            # in this case, the three quadratic subfields are all ramified, and which
            # one survives the quotient depends on the *second* inertia group:
            output.extend([(3, 2 * field.slopes[0] + (field.slopes[1] - field.slopes[0]), 4),
                           (2, 2 * field.slopes[1], 2)
                          ])

        else:
            raise Exception("Not possible")
    return output


def do_d4_fields():
    output = []
    # return output
    for field in make_d4_fields():
        # Note that our representation of d4_fields is as *quartic* fields.
        # Thus, only *half* of the automorphisms (the inner ones) keep the quadratic
        # subfield fixed. Moreover, it means that the conductor is c - d.
        d = get_quad_disc(field.deg2_subfield)
        output.append((d, field.c - d, 4))
    return output


def parse_line(line, tuple_type):
    dict_version = {'c': int(line[0]),
                    'e': int(line[1]),
                    'f': int(line[2]),
                    'd': line[3],
                    'eps': line[4],
                    'poly': line[5],
                    'G': line[6],
                    'inertia': line[7],
                    'slopes': eval(line[8])}
    if len(line) == 10:
        dict_version['deg2_subfield'] = line[9]
    return tuple_type(**dict_version)


def parse_raw(raw, tuple_type):
    reader = csv.reader(StringIO(raw))
    next(reader)
    return [parse_line(line, tuple_type) for line in reader]


def compute():
    outcomes = defaultdict(int)
    for func in (do_unram_fields, do_c2_fields, do_c4_fields, do_v4_fields, do_d4_fields):
        for (d, cond, count) in func():
            assert cond - d >= 0
            outcomes[(d, cond - d)] += count
    return outcomes


if __name__ == '__main__':
    outcomes = compute()
    print(outcomes)
    print(sum(val / 2**(d+q+3) for (d, q), val in outcomes.items()))
