#!/usr/bin/env python
import os, os.path

masterEnv = Environment(
    tools=['default','gfortran'],
    F90='gfortran',
    LINK='gfortran',
    FORTRANMODDIR='#/include',
    FORTRANMODDIRPREFIX='-J',
    F90PATH=['#/include'])

Export('masterEnv')

masterEnv.SConscript(['src/linalg/SConstruct', 'src/mesh/SConstruct', 'src/fem/SConstruct','examples/SConstruct'])
