from setuptools import setup, Extension

setup(
    name="ao",
    ext_modules=[
        Extension(
            "ao",
            sources=['aomodule.cpp'],
            language="c++",                        
            # extra_compile_args = ["-O3", "-std=c++17"],
            extra_compile_args = ["/O2", "/std:c++17",], # msvc的选项，没有O3，O2已经是启用最大速度优化
            compiler_directives={'language_level' : "3"},
        )
    ]
)