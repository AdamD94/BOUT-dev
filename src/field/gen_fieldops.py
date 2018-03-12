#!/usr/bin/env python3
"""Code-generator for arithmetic operators on Field2Ds/Field3Ds

This uses the jinja template in gen_fieldops.jinja to generate code
for the arithmetic operators, and prints to stdout.

The `Field` class provides some helper functions for determining how to
pass a variable by reference or pointer, and how to name arguments in
function signatures. This allows us to push some logic into the
templates themselves.

"""


from __future__ import print_function

try:
    from builtins import object
except ImportError:
    pass

import argparse
from collections import OrderedDict
import contextlib
from copy import deepcopy as copy
import itertools
import sys

try:
    import jinja2
except ImportError:
    raise ImportError('Missing Python module "jinja2". See "Field2D/Field3D Arithmetic '
                      'Operators" in the BOUT++ user manual for more information')


# This allows us to open either a file or stdout with the same code
@contextlib.contextmanager
def smart_open(filename, mode='r'):
    """Open stdin or stdout using a contextmanager

    From: http://stackoverflow.com/a/29824059/2043465
    """
    if filename == '-':
        if mode is None or mode == '' or 'r' in mode:
            fh = sys.stdin
        else:
            fh = sys.stdout
    else:
        fh = open(filename, mode)
    try:
        yield fh
    finally:
        if filename is not '-':
            fh.close()


# The arthimetic operators
# OrderedDict to (try to) ensure consistency between python 2 & 3
operators = OrderedDict([
    ('*', 'multiplication'),
    ('/', 'division'),
    ('+', 'addition'),
    ('-', 'subtraction'),
])

header = """// This file is autogenerated - see gen_fieldops.py
#include <bout/mesh.hxx>
#include <bout/region.hxx>

#include <field2d.hxx>
#include <field3d.hxx>
#include <globals.hxx>
#include <interpolation.hxx>
"""


class Field(object):
    """Abstracts over BoutReals and Field2D/3D/Perps

    Provides some helper functions for writing function signatures and
    passing data

    """

    def __init__(self, field_type, dimensions, name=None, index_var=None,
                 jz_var="jz", mixed_base_ind = "mixed_base_ind"):
        # C++ type of the field, e.g. Field3D
        self.field_type = field_type
        # array: dimensions of the field
        self.dimensions = dimensions
        # name of this field
        self.name = name
        # Name of the indexing variable
        self.index_var = index_var
        # Name of jz variable
        self.jz_var = jz_var
        self.mixed_base_ind = mixed_base_ind
        if self.field_type == "Field3D":
            self.region_type="3D"
        elif self.field_type == "Field2D":
            self.region_type="2D"
        else:
            self.region_type="_INVALID_"

    @property
    def passByReference(self):
        """Returns "Type& name", except if field_type is BoutReal,
        in which case just returns "Type name"

        """
        return "{self.field_type}{ref} {self.name}".format(
            self=self, ref="&" if self.field_type != "BoutReal" else "")

    @property
    def index(self):
        """Returns "[{index_var}]", except if field_type is BoutReal,
        in which case just returns ""

        """
        if self.field_type == "BoutReal":
            return "{self.name}".format(self=self)
        else:
            return "{self.name}[{self.index_var}]".format(self=self)

    @property
    def mixed_index(self):
        """Returns "[{index_var} + {jz_var}]" if field_type is Field3D,
        self.index if Field2D or just returns "" for BoutReal

        """
        if self.field_type == "BoutReal":
            return "{self.name}".format(self=self)
        elif self.field_type == "Field3D":
            return "{self.name}[{self.mixed_base_ind} + {self.jz_var}]".format(self=self)
        else:  # Field2D
            return "{self.name}[{self.index_var}]".format(self=self)

    def __eq__(self, other):
        try:
            return self.field_type == other.field_type
        except AttributeError:
            return self.field_type == other

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        return "Field({}, {}, {})".format(self.field_type, self.dimensions, self.name)

    def __str__(self):
        return self.field_type


def returnType(f1, f2):
    """Determine a suitable return type, by seeing which field is 'larger'.

    """
    if f1 == f2:
        return copy(f1)
    elif f1 == 'BoutReal':
        return copy(f2)
    elif f2 == 'BoutReal':
        return copy(f1)
    else:
        return copy(field3D)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Generate code for the Field arithmetic operators")
    # By default write to stdout
    parser.add_argument("--filename", default="-",
                        help="Write output to FILENAME instead of stdout")
    # By default use OpenMP enabled loops but allow to disable
    parser.add_argument("--no-openmp", action="store_false", default=False, dest = "noOpenMP", 
                        help="Don't use OpenMP compatible loops")

    args = parser.parse_args()

    #Setup
    index_var = 'index'
    jz_var = 'jz'
    region_name = '"RGN_ALL"'
    
    if args.noOpenMP:
        region_loop = 'BLOCK_REGION_LOOP_SERIAL'
    else:
        region_loop = 'BLOCK_REGION_LOOP'
        
    # Declare what fields we currently support:
    # Field perp is currently missing
    field3D = Field('Field3D', ['x', 'y', 'z'], index_var=index_var)
    field2D = Field('Field2D', ['x', 'y'], index_var=index_var)
    boutreal = Field('BoutReal', [], index_var=index_var)
    fields = [field3D, field2D, boutreal]

    with smart_open(args.filename, "w") as f:
        f.write(header)
        f.write("\n")

    env = jinja2.Environment(loader=jinja2.FileSystemLoader('.'),
                             trim_blocks=True)

    template = env.get_template("gen_fieldops.jinja")

    for lhs, rhs in itertools.product(fields, fields):
        # We don't have to define BoutReal BoutReal operations
        if lhs == rhs == 'BoutReal':
            continue
        rhs = copy(rhs)
        lhs = copy(lhs)

        # The output of the operation. The `larger` of the two fields.
        out = returnType(rhs, lhs)         
        out.name = 'result'
        lhs.name = 'lhs'
        rhs.name = 'rhs'

        mixed = 'no'
        region_type = out.region_type

        # Determine if we have 3D/2D or 2D/3D
        if out.field_type == "Field3D":
            if lhs.field_type == "Field2D":
                mixed = '2d3d'
                region_type = lhs.region_type
            elif rhs.field_type == "Field2D":
                mixed = '3d2d'
                region_type = rhs.region_type
                
        for operator, operator_name in operators.items():

            template_args = {
                'operator': operator,
                'operator_name': operator_name,
                
                'out': out,
                'lhs': lhs,
                'rhs': rhs,

                'mixed': mixed,
                'mixed_base_ind': 'mixed_base_ind',
                
                'region_loop': region_loop,
                'region_name': region_name,
                'region_type': region_type,
                'index_var': index_var,
                'jz_var': jz_var
                
            }

            with smart_open(args.filename, "a") as f:
                f.write(template.render(**template_args))
                f.write("\n")
