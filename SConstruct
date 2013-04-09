#!/usr/bin/env python

env = Environment(
    tools=['default','gfortran'],
    F90='gfortran',
    LINK='gfortran',
    FORTRANMODDIR='include',
    FORTRANMODDIRPREFIX='-J',
    F90PATH='include')

env.SConscript(['src/linalg/SConstruct', 'src/mesh/SConstruct', 'src/fem/SConstruct','examples/SConstruct'])
