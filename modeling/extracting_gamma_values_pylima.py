# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 20:07:11 2019

@author: rstreet
"""

>>> from pyLIMA import stars
>>> ss = stars.Star()
>>> ss.T_eff = 4250 
>>> ss.log_g = 4.5 
>>> ss.find_gamma("g'")
0.7992908463189327
>>> ss.find_gamma("r'")
0.6676151520540636
>>> ss.find_gamma("ip'")
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/ebachelet/.local/lib/python3.6/site-packages/pyLIMA-0.1.5-py3.6.egg/pyLIMA/stars.py", line 97, in find_gamma
KeyError: "ip'"
>>> ss.find_gamma("i'")
0.5418343571277271
>>> ss.find_gamma("Z")
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/ebachelet/.local/lib/python3.6/site-packages/pyLIMA-0.1.5-py3.6.egg/pyLIMA/stars.py", line 97, in find_gamma
KeyError: 'Z'
>>> ss.claret_table[:,-1]
array(['Kp', 'Kp', 'Kp', ..., "z'", "z'", "z'"], dtype='<U32')
>>> ss.find_gamma("Zp")
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/home/ebachelet/.local/lib/python3.6/site-packages/pyLIMA-0.1.5-py3.6.egg/pyLIMA/stars.py", line 97, in find_gamma
KeyError: 'Zp'
>>> np.unique(ss.claret_table[:,-1])
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
NameError: name 'np' is not defined
>>> import numpy as np
>>> np.unique(ss.claret_table[:,-1])
array(['B', 'C', 'H', 'I', 'J', 'K', 'Kp', 'R', 'S1', 'S2', 'S3', 'S4',
       'U', 'V', 'b', "g'", "i'", "r'", 'u', "u'", 'v', 'y', "z'"],
      dtype='<U32')
>>> ss.find_gamma("z'")
0.45559466317426534


