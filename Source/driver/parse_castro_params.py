#!/usr/bin/env python3

"""
This script parses the list of C++ runtime parameters and writes the
necessary header files and Fortran routines to make them available
in Castro's C++ routines and (optionally) the Fortran routines
through meth_params_module.

parameters have the format:

  name  type  default  need-in-fortran?  ifdef fortran-name  fortran-type

the first three (name, type, default) are mandatory:

  name: the name of the parameter.  This will be the same name as the
    variable in C++ unless a pair is specified as (name, cpp_name)

  type: the C++ data type (int, bool, Real, string)

  default: the default value.  If specified as a pair, (a, b), then
    the first value is the normal default and the second is for
    debug mode (#ifdef AMREX_DEBUG)

the next are optional:

   need-in-fortran: if "y" then we do a pp.query() in meth_params.F90

   ifdef: only define this parameter if the name provided is #ifdef-ed

   fortran-name: if a different variable name in Fortran, specify here

   fortran-type: if a different data type in Fortran, specify here

Any line beginning with a "#" is ignored

Commands begin with a "@":

   @namespace: sets the namespace that these will be under (see below)
     it also gives the C++ class name.
     if we include the keyword "static" after the name, then the parameters
     will be defined as static member variables in C++

     e.g. @namespace castro Castro static

Note: categories listed in the input file aren't used for code generation
but are used for the documentation generation


For a namespace, name, we write out:

  -- name_params.H  (for castro, included in Castro.H):
     declares the static variables of the Castro class

  -- name_defaults.H  (for castro, included in Castro.cpp):
     sets the defaults of the runtime parameters

  -- name_queries.H  (for castro, included in Castro.cpp):
     does the parmparse query to override the default in C++

we write out a single copy of:

  -- meth_params.F90
     does the parmparse query to override the default in Fortran,
     and sets a number of other parameters specific to the F90 routinse

"""

import argparse
import re
import sys

FWARNING = """
! This file is automatically created by parse_castro_params.py.  To update
! or add runtime parameters, please edit _cpp_parameters and then run
! mk_params.sh\n
"""

CWARNING = """
// This file is automatically created by parse_castro_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh\n
"""

param_include_dir = "param_includes/"


class Param:
    """ the basic parameter class.  For each parameter, we hold the name,
        type, and default.  For some parameters, we also take a second
        value of the default, for use in debug mode (delimited via
        #ifdef AMREX_DEBUG)

    """

    def __init__(self, name, dtype, default,
                 cpp_var_name=None,
                 namespace=None, cpp_class=None, static=None,
                 debug_default=None,
                 in_fortran=0, f90_name=None, f90_dtype=None,
                 ifdef=None):

        self.name = name
        self.dtype = dtype
        self.default = default
        self.cpp_var_name = cpp_var_name

        self.namespace = namespace
        self.cpp_class = cpp_class

        if static is None:
            self.static = 0
        else:
            self.static = static

        self.debug_default = debug_default
        self.in_fortran = in_fortran

        if ifdef == "None":
            self.ifdef = None
        else:
            self.ifdef = ifdef

        if f90_name is None:
            self.f90_name = name
        else:
            self.f90_name = f90_name

        if f90_dtype is None:
            self.f90_dtype = dtype
        else:
            self.f90_dtype = f90_dtype

    def get_default_string(self):
        # this is the line that goes into castro_defaults.H included
        # into Castro.cpp

        if self.dtype == "int":
            tstr = "int         {}::{}".format(self.cpp_class, self.cpp_var_name)
        elif self.dtype == "bool":
            tstr = "bool        {}::{}".format(self.cpp_class, self.cpp_var_name)
        elif self.dtype == "Real":
            tstr = "amrex::Real {}::{}".format(self.cpp_class, self.cpp_var_name)
        elif self.dtype == "string":
            tstr = "std::string {}::{}".format(self.cpp_class, self.cpp_var_name)
        else:
            sys.exit("invalid data type for parameter {}".format(self.name))

        ostr = ""

        if not self.debug_default is None:
            ostr += "#ifdef AMREX_DEBUG\n"
            ostr += "{} = {};\n".format(tstr, self.debug_default)
            ostr += "#else\n"
            ostr += "{} = {};\n".format(tstr, self.default)
            ostr += "#endif\n"
        else:
            ostr += "{} = {};\n".format(tstr, self.default)

        return ostr

    def get_f90_default_string(self):
        # this is the line that goes into set_castro_method_params()
        # to set the default value of the variable

        ostr = ""

        # convert to the double precision notation Fortran knows
        # if the parameter is already of the form "#.e###" then
        # it is easy as swapping out "e" for "d"; if it is a number
        # like 0.1 without a format specifier, then add a d0 to it
        # because the C++ will read it in that way and we want to
        # give identical results (at least to within roundoff)

        if self.debug_default is not None:
            debug_default = self.debug_default
            if self.dtype == "Real":
                if "d" in debug_default:
                    debug_default = debug_default.replace("d", "e")
                debug_default += "_rt"

        default = self.default
        if self.dtype == "Real":
            if "d" in default:
                default = default.replace("d", "e")
            default += "_rt"

        name = self.f90_name

        # for a character, we need to allocate its length.  We allocate
        # to 1, and the Fortran parmparse will resize
        if self.dtype == "string":
            ostr += "    allocate(character(len=1)::{})\n".format(name)
        else:
            ostr += "    allocate({})\n".format(name)

        if not self.debug_default is None:
            ostr += "#ifdef AMREX_DEBUG\n"
            ostr += "    {} = {};\n".format(name, debug_default)
            ostr += "#else\n"
            ostr += "    {} = {};\n".format(name, default)
            ostr += "#endif\n"
        else:
            ostr += "    {} = {};\n".format(name, default)

        return ostr

    def get_cuda_managed_string(self):
        """this is the string that sets the variable as managed for CUDA"""
        if self.f90_dtype == "string":
            return "\n"

        cstr = ""
        if self.ifdef is not None:
            cstr += "#ifdef {}\n".format(self.ifdef)
        cstr += "attributes(managed) :: {}\n".format(self.f90_name)
        if self.ifdef is not None:
            cstr += "#endif\n"
        return cstr

    def get_query_string(self, language):
        # this is the line that queries the ParmParse object to get
        # the value of the runtime parameter from the inputs file.
        # This goes into castro_queries.H included into Castro.cpp

        ostr = ""
        if language == "C++":
            ostr += "pp.query(\"{}\", {});\n".format(self.name, self.cpp_var_name)
        elif language == "F90":
            ostr += "    call pp%query(\"{}\", {})\n".format(self.name, self.f90_name)
        else:
            sys.exit("invalid language choice in get_query_string")

        return ostr

    def default_format(self):
        """return the variable in a format that it can be recognized in C++ code"""
        if self.dtype == "string":
            return '{}'.format(self.default)

        return self.default

    def get_job_info_test(self):
        # this is the output in C++ in the job_info writing

        ostr = 'jobInfoFile << ({}::{} == {} ? "    " : "[*] ") << "{}.{} = " << {}::{} << std::endl;\n'.format(
            self.cpp_class, self.cpp_var_name, self.default_format(),
            self.namespace, self.cpp_var_name,
            self.cpp_class, self.cpp_var_name)

        return ostr


    def get_decl_string(self):
        # this is the line that goes into castro_params.H included
        # into Castro.H

        static = ""
        if self.static:
            static = "static"

        if self.dtype == "int":
            tstr = "{} int {};\n".format(static, self.cpp_var_name)
        elif self.dtype == "bool":
            tstr = "{} bool {};\n".format(static, self.cpp_var_name)
        elif self.dtype == "Real":
            tstr = "{} amrex::Real {};\n".format(static, self.cpp_var_name)
        elif self.dtype == "string":
            tstr = "{} std::string {};\n".format(static, self.cpp_var_name)
        else:
            sys.exit("invalid data type for parameter {}".format(self.name))

        ostr = ""
        ostr += tstr

        return ostr

    def get_f90_decl_string(self):
        # this is the line that goes into meth_params.f90

        if not self.in_fortran:
            return None

        if self.f90_dtype == "int":
            tstr = "integer,  allocatable, save :: {}\n".format(self.f90_name)
        elif self.f90_dtype == "Real":
            tstr = "real(rt), allocatable, save :: {}\n".format(self.f90_name)
        elif self.f90_dtype == "logical":
            tstr = "logical,  allocatable, save :: {}\n".format(self.f90_name)
        elif self.f90_dtype == "string":
            tstr = "character (len=:), allocatable, save :: {}\n".format(self.f90_name)
            print("warning: string parameter {} will not be available on the GPU".format(
                self.f90_name))
        else:
            sys.exit("unsupported datatype for Fortran: {}".format(self.name))

        return tstr


def write_meth_module(plist, meth_template):
    """this writes the meth_params_module, starting with the meth_template
       and inserting the runtime parameter declaration in the correct
       place
    """

    try:
        mt = open(meth_template, "r")
    except IOError:
        sys.exit("invalid template file")

    try:
        mo = open("meth_params.F90", "w")
    except IOError:
        sys.exit("unable to open meth_params.F90 for writing")


    mo.write(FWARNING)

    param_decls = [p.get_f90_decl_string() for p in plist if p.in_fortran == 1]
    params = [p for p in plist if p.in_fortran == 1]

    decls = ""

    for p in param_decls:
        decls += "  {}".format(p)

    cuda_managed_decls = [p.get_cuda_managed_string() for p in plist if p.in_fortran == 1]

    cuda_managed_string = ""
    for p in cuda_managed_decls:
        cuda_managed_string += "{}".format(p)

    for line in mt:
        if line.find("@@f90_declarations@@") > 0:
            mo.write(decls)

            # Do the CUDA managed declarations

            mo.write("\n")
            mo.write("#ifdef AMREX_USE_CUDA\n")
            mo.write(cuda_managed_string)
            mo.write("#endif\n")

            # Now do the OpenACC declarations

            mo.write("\n")
            mo.write("  !$acc declare &\n")
            for n, p in enumerate(params):
                if p.f90_dtype == "string":
                    print("warning: string parameter {} will not be on the GPU".format(p.name),
                          file=sys.stderr)
                    continue

                if p.ifdef is not None:
                    mo.write("#ifdef {}\n".format(p.ifdef))
                mo.write("  !$acc create({})".format(p.f90_name))

                if n != len(params)-1:
                    mo.write(" &\n")
                else:
                    mo.write("\n")

                if p.ifdef is not None:
                    mo.write("#endif\n")



        elif line.find("@@set_castro_params@@") >= 0:

            namespaces = {q.namespace for q in params}
            print("Fortran namespaces: ", namespaces)
            for nm in namespaces:
                params_nm = [q for q in params if q.namespace == nm]
                ifdefs = {q.ifdef for q in params_nm}

                for ifdef in ifdefs:
                    if ifdef is None:
                        for p in [q for q in params_nm if q.ifdef is None]:
                            mo.write(p.get_f90_default_string())
                    else:
                        mo.write("#ifdef {}\n".format(ifdef))
                        for p in [q for q in params_nm if q.ifdef == ifdef]:
                            mo.write(p.get_f90_default_string())
                        mo.write("#endif\n")

                mo.write("\n")

                mo.write('    call amrex_parmparse_build(pp, "{}")\n'.format(nm))

                for ifdef in ifdefs:
                    if ifdef is None:
                        for p in [q for q in params_nm if q.ifdef is None]:
                            mo.write(p.get_query_string("F90"))
                    else:
                        mo.write("#ifdef {}\n".format(ifdef))
                        for p in [q for q in params_nm if q.ifdef == ifdef]:
                            mo.write(p.get_query_string("F90"))
                        mo.write("#endif\n")

                mo.write('    call amrex_parmparse_destroy(pp)\n')

                mo.write("\n\n")

            # Now do the OpenACC device updates

            mo.write("\n")
            mo.write("    !$acc update &\n")
            for n, p in enumerate(params):
                if p.f90_dtype == "string": continue

                if p.ifdef is not None:
                    mo.write("#ifdef {}\n".format(p.ifdef))
                mo.write("    !$acc device({})".format(p.f90_name))

                if n != len(params)-1:
                    mo.write(" &\n")
                else:
                    mo.write("\n")

                if p.ifdef is not None:
                    mo.write("#endif\n")

        elif line.find("@@free_castro_params@@") >= 0:

            params_free = [q for q in params if q.in_fortran == 1]

            for p in params_free:
                mo.write("    if (allocated({})) then\n".format(p.f90_name))
                mo.write("        deallocate({})\n".format(p.f90_name))
                mo.write("    end if\n")

            mo.write("\n\n")


        else:
            mo.write(line)

    mo.close()
    mt.close()


def parse_params(infile, meth_template):

    params = []

    namespace = None
    cpp_class = None
    static = None

    try:
        f = open(infile)
    except IOError:
        sys.exit("error openning the input file")


    for line in f:
        if line[0] == "#":
            continue

        if line.strip() == "":
            continue

        if line[0] == "@":
            # this is a command
            cmd, value = line.split(":")
            if cmd == "@namespace":
                fields = value.split()
                namespace = fields[0]
                cpp_class = fields[1]

                try:
                    static = fields[2]
                except IndexError:
                    static = ""

                # do we have the static keyword?
                if "static" in static:
                    static = 1
                else:
                    static = 0

            else:
                sys.exit("invalid command")

            continue

        # this splits the line into separate fields.  A field is a
        # single word or a pair in parentheses like "(a, b)"
        fields = re.findall(r'[\w\"\+\.-]+|\([\w+\.-]+\s*,\s*[\w\+\.-]+\)', line)

        name = fields[0]
        if name[0] == "(":
            name, cpp_var_name = re.findall(r"\w+", name)
        else:
            cpp_var_name = name

        dtype = fields[1]

        default = fields[2]
        if default[0] == "(":
            default, debug_default = re.findall(r"\w+", default)
        else:
            debug_default = None

        try:
            in_fortran_string = fields[3]
        except IndexError:
            in_fortran = 0
        else:
            if in_fortran_string.lower().strip() == "y":
                in_fortran = 1
            else:
                in_fortran = 0

        try:
            ifdef = fields[4]
        except IndexError:
            ifdef = None

        try:
            f90_name = fields[5]
        except IndexError:
            f90_name = None

        try:
            f90_dtype = fields[6]
        except IndexError:
            f90_dtype = None

        if namespace is None:
            sys.exit("namespace not set")

        params.append(Param(name, dtype, default,
                            cpp_var_name=cpp_var_name,
                            namespace=namespace,
                            cpp_class=cpp_class,
                            static=static,
                            debug_default=debug_default,
                            in_fortran=in_fortran, f90_name=f90_name, f90_dtype=f90_dtype,
                            ifdef=ifdef))



    # output

    # find all the namespaces
    namespaces = {q.namespace for q in params}

    for nm in namespaces:

        params_nm = [q for q in params if q.namespace == nm]
        ifdefs = {q.ifdef for q in params_nm}

        # write name_defaults.H
        try:
            cd = open("{}/{}_defaults.H".format(param_include_dir, nm), "w")
        except IOError:
            sys.exit("unable to open {}_defaults.H for writing".format(nm))

        cd.write(CWARNING)

        for ifdef in ifdefs:
            if ifdef is None:
                for p in [q for q in params_nm if q.ifdef is None]:
                    cd.write(p.get_default_string())
            else:
                cd.write("#ifdef {}\n".format(ifdef))
                for p in [q for q in params_nm if q.ifdef == ifdef]:
                    cd.write(p.get_default_string())
                cd.write("#endif\n")

        cd.close()

        # write name_params.H
        try:
            cp = open("{}/{}_params.H".format(param_include_dir, nm), "w")
        except IOError:
            sys.exit("unable to open {}_params.H for writing".format(nm))

        cp.write(CWARNING)

        for ifdef in ifdefs:
            if ifdef is None:
                for p in [q for q in params_nm if q.ifdef is None]:
                    cp.write(p.get_decl_string())
            else:
                cp.write("#ifdef {}\n".format(ifdef))
                for p in [q for q in params_nm if q.ifdef == ifdef]:
                    cp.write(p.get_decl_string())
                cp.write("#endif\n")

        cp.close()

        # write castro_queries.H
        try:
            cq = open("{}/{}_queries.H".format(param_include_dir, nm), "w")
        except IOError:
            sys.exit("unable to open {}_queries.H for writing".format(nm))

        cq.write(CWARNING)

        for ifdef in ifdefs:
            if ifdef is None:
                for p in [q for q in params_nm if q.ifdef is None]:
                    cq.write(p.get_query_string("C++"))
            else:
                cq.write("#ifdef {}\n".format(ifdef))
                for p in [q for q in params_nm if q.ifdef == ifdef]:
                    cq.write(p.get_query_string("C++"))
                cq.write("#endif\n")

        cq.close()

        # write the job info tests
        try:
            jo = open("{}/{}_job_info_tests.H".format(param_include_dir, nm), "w")
        except IOError:
            sys.exit("unable to open {}_job_info_tests.H".format(nm))

        for ifdef in ifdefs:
            if ifdef is None:
                for p in [q for q in params_nm if q.ifdef is None]:
                    jo.write(p.get_job_info_test())
            else:
                jo.write("#ifdef {}\n".format(ifdef))
                for p in [q for q in params_nm if q.ifdef == ifdef]:
                    jo.write(p.get_job_info_test())
                jo.write("#endif\n")

        jo.close()

    # write the Fortran module
    write_meth_module(params, meth_template)


def main():
    """the main driver"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", type=str, default=None,
                        help="template for the meth_params module")
    parser.add_argument("input_file", type=str, nargs=1,
                        help="input file containing the list of parameters we will define")

    args = parser.parse_args()

    parse_params(args.input_file[0], args.m)

if __name__ == "__main__":
    main()
