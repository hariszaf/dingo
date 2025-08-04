# dingo : a python library for metabolic networks sampling and analysis
# dingo is part of GeomScale project

# Copyright (c) 2021 Apostolos Chalkis
# Copyright (c) 2024 Vissarion Fisikopoulos

# Licensed under GNU LGPL.3, see LICENCE file

# This is the setup Python script for building the dingo library
import os
import numpy
import platform

from os.path import join
from Cython.Build import cythonize
from setuptools import setup, Extension

# Determine the OS
current_platform = platform.system()

# Compiler arguments
# link_args              = ["-O3", "-fopenmp"]
# compiler_args          = ["-std=c++17", "-O3", "-DBOOST_NO_AUTO_PTR", "-ldl", "-lm", "-fopenmp"]
# lp_solve_compiler_args = ["-DYY_NEVER_INTERACTIVE", "-DLoadInverseLib=0", "-DLoadLanguageLib=0",
#                           "-DRoleIsExternalInvEngine", "-DINVERSE_ACTIVE=3", "-DLoadableBlasLib=0"]

# =========================
# # Add these flags to explicitly disable SIMD extensions on x86
# disable_simd_flags = ["-mno-sse", "-mno-sse2", "-mno-avx"]

# # Compiler arguments
# link_args = ["-O3"]
# compiler_args = ["-std=c++17", "-O3", "-DBOOST_NO_AUTO_PTR", "-ldl", "-lm"]
# lp_solve_compiler_args = [
#     "-DYY_NEVER_INTERACTIVE", "-DLoadInverseLib=0", "-DLoadLanguageLib=0",
#     "-DRoleIsExternalInvEngine", "-DINVERSE_ACTIVE=3", "-DLoadableBlasLib=0"
# ]

# # Set specific arguments for Linux or macOS
# if current_platform == "Linux":
#     link_args.append("-fopenmp")
#     compiler_args.append("-fopenmp")
# elif current_platform == "Darwin":  # macOS
#     link_args.extend(["-Xpreprocessor", "-fopenmp", "-lomp"])
#     compiler_args.extend(["-Xpreprocessor", "-fopenmp", "-lomp"])

# # Apply SIMD-disabling flags only on x86 systems
# arch = platform.machine()
# if arch in ("x86_64", "i386", "i686"):
#     compiler_args = disable_simd_flags + compiler_args
# ===========================

# Determine platform
current_platform = platform.system()
arch = platform.machine()

# Start with generic settings
base_link_args = ["-O3", "-fopenmp", "-ldl", "-lm"]
lp_solve_compiler_args = [
    "-DYY_NEVER_INTERACTIVE",
    "-DLoadInverseLib=0",
    "-DLoadLanguageLib=0",
    "-DRoleIsExternalInvEngine",
    "-DINVERSE_ACTIVE=3",
    "-DLoadableBlasLib=0"
]
base_compiler_args = ["-std=c++17", "-O3", "-fopenmp"]

# SIMD disabling (only for x86 architectures)
disable_simd_flags = []

# Platform-specific settings
# std::is_same_v --> c++17

if current_platform == "Darwin":
    # if arch in ("x86_64", "i386", "i686"):
    #     disable_simd_flags = ["-mno-sse", "-mno-sse2", "-mno-avx"]

    # Get OpenMP paths from brew
    brew_prefix = os.popen("brew --prefix libomp").read().strip()
    omp_include = os.path.join(brew_prefix, "include")
    omp_lib = os.path.join(brew_prefix, "lib")

    base_compiler_args.extend([
        "-Xpreprocessor",
        "-stdlib=libc++",
        f"-I{omp_include}"
    ])
    base_link_args.extend([
        "-Xpreprocessor",
        "-lomp",
        "-stdlib=libc++",
        f"-L{omp_lib}"
    ])
else:
    base_compiler_args.append("-DBOOST_NO_AUTO_PTR")

# Final flags
compiler_args = disable_simd_flags + base_compiler_args + lp_solve_compiler_args
link_args = base_link_args

# Ext
volesti_include_dirs = [
    # include binding files
    join("dingo", "bindings"),
    # the volesti code uses some external classes.
    # external directories we need to add
    join("eigen"),
    join("boost_1_76_0"),
    join("boost_1_76_0", "boost"),
    join("lp_solve_5.5"),
    join("lp_solve_5.5", "bfp"),
    join("lp_solve_5.5", "bfp", "bfp_LUSOL"),
    join("lp_solve_5.5", "bfp", "bfp_LUSOL", "LUSOL"),
    join("lp_solve_5.5", "colamd"),
    join("lp_solve_5.5", "shared"),
    join("volesti", "external"),
    join("volesti", "external", "minimum_ellipsoid"),
    # include and add the directories on the "include" directory
    join("volesti", "include"),
    join("volesti", "include", "convex_bodies"),
    join("volesti", "include", "random_walks"),
    join("volesti", "include", "volume"),
    join("volesti", "include", "generators"),
    join("volesti", "include", "cartesian_geom"),
]

src_files = [
    "lp_solve_5.5/bfp/bfp_LUSOL/lp_LUSOL.c",
    "lp_solve_5.5/bfp/bfp_LUSOL/LUSOL/lusol.c",
    "lp_solve_5.5/colamd/colamd.c",
    "lp_solve_5.5/ini.c",
    "lp_solve_5.5/shared/commonlib.c",
    "lp_solve_5.5/shared/mmio.c",
    "lp_solve_5.5/shared/myblas.c",
    "lp_solve_5.5/lp_crash.c",
    "lp_solve_5.5/lp_Hash.c",
    "lp_solve_5.5/lp_lib.c",
    "lp_solve_5.5/lp_matrix.c",
    "lp_solve_5.5/lp_MDO.c",
    "lp_solve_5.5/lp_mipbb.c",
    "lp_solve_5.5/lp_MPS.c",
    "lp_solve_5.5/lp_params.c",
    "lp_solve_5.5/lp_presolve.c",
    "lp_solve_5.5/lp_price.c",
    "lp_solve_5.5/lp_pricePSE.c",
    "lp_solve_5.5/lp_report.c",
    "lp_solve_5.5/lp_scale.c",
    "lp_solve_5.5/lp_simplex.c",
    "lp_solve_5.5/lp_SOS.c",
    "lp_solve_5.5/lp_utils.c",
    "lp_solve_5.5/lp_wlp.c",
    "dingo/volestipy.pyx",
    "dingo/bindings/bindings.cpp"
]

# Return the directory that contains the NumPy *.h header files.
# Extension modules that need to compile against NumPy should use this
# function to locate the appropriate include directory.
numpy_dirs       = [numpy.get_include()]
suitesparse_dirs = ["/usr/include/suitesparse"]  # Include the SuiteSparse headers
include_dirs     = volesti_include_dirs + suitesparse_dirs + numpy_dirs

# --- Extension ---
print("Using compiler args:", compiler_args)
print("Using linker args:", link_args)

volesti_module = Extension(
    name               = "dingo.volestipy",
    language           = "c++",
    sources            = src_files,
    include_dirs       = include_dirs,
    extra_compile_args = compiler_args,
    extra_link_args    = link_args,
)

ext_modules = cythonize(
    [volesti_module],
    gdb_debug=False
)


if __name__ == "__main__":
    setup(
        packages     = ["dingo", "dingo.bindings"],
        ext_modules  = ext_modules,
        zip_safe     = False,
    )
