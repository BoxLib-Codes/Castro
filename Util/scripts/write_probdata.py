#!/usr/bin/env python3

"""This routine parses plain-text parameter files that list runtime
parameters for use in our codes.  The general format of a parameter
is:

max_step                            integer            1
small_dt                            real               1.d-10
xlo_boundary_type                   character          ""
octant                              logical            .false.

This specifies the parameter name, datatype, and default
value.

The optional 4th column indicates whether the parameter appears
in the fortin namelist ("y" or "Y").

If the parameter is an array, the optional 5th column indicates the size of
the array.

This script takes a template file and replaces keywords in it
(delimited by @@...@@) with the Fortran code required to
initialize the parameters, setup a namelist, set the defaults, etc.

Note: there are two types of parameters here, the ones that are in the
namelist are true runtime parameters.  The non-namelist parameters
should be avoided if at all possible.

"""

import argparse
import os
import re
import sys

import runtime_parameters as rp

HEADER = """
! DO NOT EDIT THIS FILE!!!
!
! This file is automatically generated by write_probdata.py at
! compile-time.
!
! To add a runtime parameter, do so by editting the appropriate _prob_params
! file.

"""

CXX_F_HEADER = """
#ifndef problem_parameters_F_H
#define problem_parameters_F_H
#include <AMReX.H>
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

void probdata_init(const int* name, const int* namlen);

void prob_params_pretty_print(int* jobinfo_file_name, const int* jobinfo_file_length);

void update_prob_params_after_cxx();

"""

CXX_F_FOOTER = """
#ifdef __cplusplus
}
#endif

#endif
"""

CXX_HEADER = """
#ifndef problem_parameters_H
#define problem_parameters_H
#include <AMReX_BLFort.H>

#include <network_properties.H>

"""

CXX_FOOTER = """
#endif
"""


def get_next_line(fin):
    """return the next, non-blank line, with comments stripped"""
    line = fin.readline()

    pos = line.find("#")

    while (pos == 0 or line.strip() == "") and line:
        line = fin.readline()
        pos = line.find("#")

    return line[:pos]


def parse_param_file(params_list, param_file):
    """read all the parameters in the prob_param_files and add valid
    parameters to the params list.  This returns the parameter list.

    """

    namespace = "problem"

    try:
        f = open(param_file, "r")
    except FileNotFoundError:
        sys.exit(f"write_probdata.py: ERROR: file {param_file} does not exist")

    line = get_next_line(f)

    err = 0

    while line and not err:

        # this splits the line into separate fields.  A field is a
        # single word or a pair in parentheses like "(a, b)"
        fields = re.findall(r'[\w\"\+\./\-]+|\([\w+\./\-]+\s*,\s*[\w\+\.\-]+\)', line)

        if len(fields) < 3:
            print("write_probdata.py: ERROR: missing one or more fields in parameter definition.")
            err = 1
            continue

        name = fields[0]
        dtype = fields[1]
        default = fields[2]

        current_param = rp.Param(name, dtype, default,
                                 namespace=namespace,
                                 in_fortran=1)

        # optional field: in namelist
        try:
            in_namelist_in = fields[3]
            if in_namelist_in in ["y", "Y"]:
                in_namelist = True
            else:
                in_namelist = False

        except IndexError:
            in_namelist = False


        # optional field: size
        try:
            size = fields[4]
        except IndexError:
            size = 1

        current_param.in_namelist = in_namelist
        current_param.size = size

        # check to see if this parameter is defined in the current
        # list if we delete the old one and take the new one (we
        # assume that later files automatically have higher
        # priority)
        p_names = [p.name for p in params_list]
        try:
            idx = p_names.index(current_param.name)
        except ValueError:
            pass
        else:
            params_list.pop(idx)

        if not err == 1:
            params_list.append(current_param)

        line = get_next_line(f)

    return err


def abort(outfile):
    """ abort exits when there is an error.  A dummy stub file is written
    out, which will cause a compilation failure """

    fout = open(outfile, "w")
    fout.write("There was an error parsing the parameter files")
    fout.close()
    sys.exit(1)


def write_probin(probin_template, prob_param_files,
                 out_file, cxx_prefix):
    """write_probin will read through the list of parameter files and
    output the new out_file"""

    params = []

    print(" ")
    print(f"write_probdata.py: creating {out_file}")

    # read the parameters defined in the parameter files

    for f in prob_param_files:
        err = parse_param_file(params, f)
        if err:
            abort(out_file)

    # open up the template
    try:
        ftemplate = open(probin_template, "r")
    except IOError:
        sys.exit(f"write_probdata.py: ERROR: file {probin_template} does not exist")

    template_lines = ftemplate.readlines()

    ftemplate.close()

    # output the template, inserting the parameter info in between the @@...@@
    fout = open(out_file, "w")

    fout.write(HEADER)

    exists_namelist = any([q.in_namelist for q in params])

    for line in template_lines:

        index = line.find("@@")

        if index >= 0:
            index2 = line.rfind("@@")

            keyword = line[index+len("@@"):index2]
            indent = index*" "

            if keyword == "declarations":

                if not exists_namelist:
                    # we always make sure there is atleast one variable in the namelist
                    fout.write(f"{indent}integer, save, public :: a_pp_dummy_var = 0\n")

                # declaraction statements for both namelist and non-namelist variables
                for p in params:
                    fout.write(f"{indent}{p.get_f90_decl_string()}")

            elif keyword == "allocations":
                for p in params:
                    fout.write(p.get_f90_default_string())

            elif keyword == "namelist":
                if not exists_namelist:
                    fout.write(f"{indent}namelist /fortin/ a_pp_dummy_var\n")
                else:
                    for p in params:
                        if p.in_namelist:
                            fout.write(f"{indent}namelist /fortin/ {p.name}\n")

            elif keyword == "printing":

                fout.write("100 format (1x, a3, 2x, a32, 1x, \"=\", 1x, a)\n")
                fout.write("101 format (1x, a3, 2x, a32, 1x, \"=\", 1x, i10)\n")
                fout.write("102 format (1x, a3, 2x, a32, 1x, \"=\", 1x, g20.10)\n")
                fout.write("103 format (1x, a3, 2x, a32, 1x, \"=\", 1x, l)\n")

                for p in params:
                    if not p.in_namelist:
                        continue

                    if p.dtype in ["bool", "logical"]:
                        ltest = f"\n{indent}ltest = {p.name} .eqv. {p.default}\n"
                    else:
                        ltest = f"\n{indent}ltest = {p.name} == {p.default}\n"

                    fout.write(ltest)

                    cmd = "merge(\"   \", \"[*]\", ltest)"

                    if p.dtype == "real":
                        fout.write(f"{indent}write (unit,102) {cmd}, &\n \"{p.name}\", {p.name}\n")

                    elif p.dtype == "string":
                        fout.write(f"{indent}write (unit,100) {cmd}, &\n \"{p.name}\", trim({p.name})\n")

                    elif p.dtype == "integer":
                        fout.write(f"{indent}write (unit,101) {cmd}, &\n \"{p.name}\", {p.name}\n")

                    elif p.dtype in ["bool", "logical"]:
                        fout.write(f"{indent}write (unit,103) {cmd}, &\n \"{p.name}\", {p.name}\n")

                    else:
                        print(f"write_probdata.py: invalid datatype for variable {p.name}")


            elif keyword == "cxx_gets":
                # this writes out the Fortran functions that can be
                # called from C++ to get the value of the parameters.

                for p in params:
                    fout.write(p.get_f90_get_function())

            elif keyword == "cxx_sets":
                # this writes out the Fortran functions that can be called from C++
                # to set the value of the parameters.

                for p in params:
                    if p.dtype in ["bool", "logical"] or p.dtype == "string":
                        continue

                    fout.write(f"  subroutine set_f90_{p.name}({p.name}_in) bind(C, name=\"set_f90_{p.name}\")\n")
                    if p.is_array():
                        fout.write(f"     {p.get_f90_decl()}, intent(in) :: {p.name}_in({p.size})\n")
                    else:
                        fout.write(f"     {p.get_f90_decl()}, intent(in) :: {p.name}_in\n")
                    fout.write(f"     {p.name} = {p.name}_in\n")
                    fout.write(f"  end subroutine set_f90_{p.name}\n\n")

            elif keyword == "fortran_parmparse_overrides":

                fout.write(f'    call amrex_parmparse_build(pp, "problem")\n')

                for p in params:
                    if p.in_namelist:
                        fout.write(p.get_query_string("F90"))

                fout.write('    call amrex_parmparse_destroy(pp)\n')
                fout.write("\n\n")

        else:
            fout.write(line)

    print(" ")
    fout.close()

    # now handle the C++ -- we need to write a header and a .cpp file
    # for the parameters + a _F.H file for the Fortran communication

    # first the _F.H file
    ofile = f"{cxx_prefix}_parameters_F.H"
    with open(ofile, "w") as fout:
        fout.write(CXX_F_HEADER)

        for p in params:
            if p.dtype == "string":
                fout.write(f"  void get_f90_{p.name}(char* {p.name});\n\n")
                fout.write(f"  void get_f90_{p.name}_len(int& slen);\n\n")
                fout.write(f"  void set_f90_{p.name}_len(int& slen);\n\n")

            else:
                fout.write(f"  void get_f90_{p.name}({p.get_cxx_decl()}* {p.name});\n\n")
                fout.write(f"  void set_f90_{p.name}({p.get_cxx_decl()}* {p.name});\n\n")

        fout.write(CXX_F_FOOTER)

    # now the main C++ header with the global data
    cxx_base = os.path.basename(cxx_prefix)

    ofile = f"{cxx_prefix}_parameters.H"
    with open(ofile, "w") as fout:
        fout.write(CXX_HEADER)

        fout.write(f"  void init_{cxx_base}_parameters();\n\n")

        fout.write(f"  void cxx_to_f90_{cxx_base}_parameters();\n\n")

        fout.write("  namespace problem {\n\n")

        for p in params:
            if p.dtype == "string":
                fout.write(f"  extern std::string {p.name};\n\n")
            else:
                if p.is_array():
                    if p.size == "nspec":
                        fout.write(f"  extern AMREX_GPU_MANAGED {p.get_cxx_decl()} {p.name}[NumSpec];\n\n")
                    else:
                        fout.write(f"  extern AMREX_GPU_MANAGED {p.get_cxx_decl()} {p.name}[{p.size}];\n\n")
                else:
                    fout.write(f"  extern AMREX_GPU_MANAGED {p.get_cxx_decl()} {p.name};\n\n")

        fout.write("  }\n\n")

        fout.write(CXX_FOOTER)

    # finally the C++ initialization routines
    ofile = f"{cxx_prefix}_parameters.cpp"
    with open(ofile, "w") as fout:
        fout.write(f"#include <{cxx_base}_parameters.H>\n")
        fout.write(f"#include <{cxx_base}_parameters_F.H>\n\n")
        fout.write("#include <AMReX_ParmParse.H>\n\n")

        for p in params:
            if p.dtype == "string":
                fout.write(f"  std::string problem::{p.name};\n\n")
            else:
                if p.is_array():
                    if p.size == "nspec":
                        fout.write(f"  AMREX_GPU_MANAGED {p.get_cxx_decl()} problem::{p.name}[NumSpec];\n\n")
                    else:
                        fout.write(f"  AMREX_GPU_MANAGED {p.get_cxx_decl()} problem::{p.name}[{p.size}];\n\n")
                else:
                    fout.write(f"  AMREX_GPU_MANAGED {p.get_cxx_decl()} problem::{p.name};\n\n")

        fout.write("\n")
        fout.write(f"  void init_{cxx_base}_parameters() {{\n")

        # first write the "get" routines to get the parameter from the
        # Fortran read -- this will either be the default or the value
        # from the probin

        fout.write("    // get the values of the parameters from Fortran\n\n")

        for p in params:
            if p.dtype == "string":
                fout.write(f"    int slen_{p.name} = 0;\n")
                fout.write(f"    get_f90_{p.name}_len(slen_{p.name});\n")
                fout.write(f"    char _{p.name}[slen_{p.name}+1];\n")
                fout.write(f"    get_f90_{p.name}(_{p.name});\n")
                fout.write(f"    problem::{p.name} = std::string(_{p.name});\n\n")
            else:
                if p.is_array():
                    fout.write(f"    get_f90_{p.name}(problem::{p.name});\n\n")
                else:
                    fout.write(f"    get_f90_{p.name}(&problem::{p.name});\n\n")


        # now write the parmparse code to get the value from the C++
        # inputs.  This will overwrite the Fortran value.

        fout.write("    // get the value from the inputs file (this overwrites the Fortran value)\n\n")

        # open namespace
        fout.write("    {\n")
        fout.write(f"      amrex::ParmParse pp(\"problem\");\n")
        for p in params:
            if p.in_namelist:
                qstr = p.get_query_string("C++")
                fout.write(f"      {qstr}")
        fout.write("    }\n")

        fout.write("  }\n")

        # finally write the set code that is called after a problem
        # does any C++ initialization which may have modified values
        # -- this sends the values back to Fortran.

        fout.write(f"  void cxx_to_f90_{cxx_base}_parameters() {{\n")
        fout.write("    int slen = 0;\n\n")

        for p in params:
            if p.dtype in ["bool", "logical"] or p.dtype == "string":
                continue

            if p.is_array():
                fout.write(f"    set_f90_{p.name}(problem::{p.name});\n\n")
            else:
                fout.write(f"    set_f90_{p.name}(&problem::{p.name});\n\n")

        fout.write("  }\n")


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', type=str, help='probin_template')
    parser.add_argument('-o', type=str, help='out_file')
    parser.add_argument('-p', type=str, help='problem parameter file names (space separated in quotes)')
    parser.add_argument('--cxx_prefix', type=str, default="prob",
                        help="a name to use in the C++ file names")

    args = parser.parse_args()

    probin_template = args.t
    out_file = args.o
    prob_params = args.p.split()

    if probin_template == "" or out_file == "":
        sys.exit("write_probdata.py: ERROR: invalid calling sequence")

    write_probin(probin_template, prob_params, out_file, args.cxx_prefix)

if __name__ == "__main__":
    main()
