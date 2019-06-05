"""
# authors: eablonet
# date : 2019
"""
import Stack as st

import matplotlib.pyplot as plt
import numpy as np
s = st.Stack()

s.read_by_path('06-04-2018', 5)
# s.select_folder()
time_nuc = int(s.datas['t_nuc_calc'].iloc[0])
print(time_nuc)
s.read_image(time_nuc-1)
s.display_contour()

s.read_image(time_nuc)
s.display_contour()
