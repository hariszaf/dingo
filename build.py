import os

# See if Cython is installed
try:
    from Cython.Build import cythonize
# Do nothing if Cython is not available
except ImportError:
    # Provide a placeholder function for the setup process if Cython is missing
    def build(setup_kwargs):
        pass
else:
    from setuptools import Extension
    from setuptools.dist import Distribution
    from distutils.command.build_ext import build_ext

    # This function will be executed in setup.py (or by setuptools directly)
    def build(setup_kwargs):
        # The file you want to compile
        extensions = ["dingo/volestipy.pyx"]

        # gcc arguments hack: enable optimizations
        os.environ["CFLAGS"] = [
            "-std=c++17",
            "-O3",
            "-DBOOST_NO_AUTO_PTR",
            "-ldl",
            "-lm",
        ]

        # Update setup kwargs to include Cython extensions and custom build
        setup_kwargs.update(
            {
                "ext_modules": cythonize(
                    extensions,
                    language_level=3,
                    compiler_directives={"linetrace": True},
                ),
                "cmdclass": {"build_ext": build_ext},
            }
        )
