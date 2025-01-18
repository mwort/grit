#!/usr/bin/env python
"""This is an abstraction module for the commandline interface used throughout this project."""

import fire  # see: https://github.com/google/python-fire/blob/master/docs/guide.md

interface = fire.Fire

if __name__ == "__main__":
    interface(__import__)
