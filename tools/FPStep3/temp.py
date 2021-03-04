#!/usr/bin/env python3
# -*-coding:Utf-8 -*
import os
from tempfile import gettempdir
from shutil import rmtree
tmp = os.path.join(gettempdir(), '.{}'.format(hash(os.times())))
os.makedirs(tmp)
print(tmp)

rmtree(tmp, ignore_errors=True)